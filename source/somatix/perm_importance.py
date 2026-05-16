import gc
import os

import numpy as np
import pandas as pd
import torch

from somatix.predict_DNN_somatic import (
    DNASeqCoverageModel_somatic,
    FEATURE_ORDER,
    annotate_variant_rows,
    convert_to_vcf,
    get_contigs_from_bam_or_default,
    load_shard_records,
    load_variant_table,
    normalize_target_chr,
    predict_feature_batch,
    shard_files,
)


FEATURE_GROUPS = {
    "base_fraction": [
        "one_hot_encoded_sequence",
        "ref_base_fractions_forward",
        "ref_base_fractions_reverse",
        "alt_base_fractions_forward",
        "alt_base_fractions_reverse",
    ],
    "target_position_coverage": ["target_position_coverage"],
    "coverage": [
        "ref_coverage_forward",
        "ref_coverage_reverse",
        "alt_coverage_forward",
        "alt_coverage_reverse",
    ],
    "mapq": [
        "ref_mapq_forward",
        "ref_mapq_reverse",
        "alt_mapq_forward",
        "alt_mapq_reverse",
    ],
    "baseq": [
        "ref_baseq_forward",
        "ref_baseq_reverse",
        "alt_baseq_forward",
        "alt_baseq_reverse",
    ],
    "mismatch_rate": [
        "ref_mismatch_rate_forward",
        "ref_mismatch_rate_reverse",
        "alt_mismatch_rate_forward",
        "alt_mismatch_rate_reverse",
    ],
}


FEATURE_GROUP_DESCRIPTIONS = {
    "base_fraction": "USED. Five-channel base-fraction branch input: padded sequence context plus ref/alt A/C/G/T/DEL/INS base-fraction pileups on forward and reverse strands; one sequence tensor shape (4, 61) and four base-fraction tensors, each shape (6, 61).",
    "target_position_coverage": "USED. Normalized total coverage at the candidate position; tensor shape per sample: (1,).",
    "coverage": "USED. Ref/alt allele coverage pileups on forward and reverse strands; four tensors, each shape per sample: (61,).",
    "mapq": "USED. Mean mapping quality for ref/alt-supporting reads on forward and reverse strands; four scalar tensors, each shape per sample: (1,).",
    "baseq": "USED. Mean base quality for ref/alt-supporting reads on forward and reverse strands; four tensors, each shape per sample: (61,).",
    "mismatch_rate": "USED. Mismatch-rate summaries for ref/alt-supporting reads on forward and reverse strands; four scalar tensors, each shape per sample: (1,).",
}


def validate_feature_groups():
    missing = {
        group_name: [name for name in feature_names if name not in FEATURE_ORDER]
        for group_name, feature_names in FEATURE_GROUPS.items()
    }
    missing = {group_name: names for group_name, names in missing.items() if names}
    if missing:
        details = "; ".join(
            f"{group_name}: {', '.join(names)}"
            for group_name, names in missing.items()
        )
        raise ValueError(f"Permutation feature group contains unsupported model features: {details}")


validate_feature_groups()


def feature_group_help_text():
    lines = [
        "Experimental feature-testing command, not the standard SomatiX prediction workflow.",
        "Available --feature groups:",
    ]
    for name, members in FEATURE_GROUPS.items():
        lines.append(f"  {name}: {FEATURE_GROUP_DESCRIPTIONS[name]}")
        lines.append(f"    datasets: {', '.join(members)}")
    return "\n".join(lines)


def check_feature_name(feature_name):
    if feature_name not in FEATURE_GROUPS:
        raise ValueError(
            f"Unknown feature group: {feature_name}\n"
            f"Available feature groups: {', '.join(FEATURE_GROUPS)}"
        )


def load_model(model_path, device):
    model = DNASeqCoverageModel_somatic().to(device)
    state = torch.load(model_path, map_location=device)
    if isinstance(state, dict) and "state_dict" in state:
        state = state["state_dict"]
    try:
        model.load_state_dict(state)
    except RuntimeError as exc:
        raise RuntimeError(
            "Failed to load model checkpoint into the current SomatiX architecture. "
            "The current model excludes total_fraction, read_length, and "
            "target_relative_position, so checkpoints trained with the old 1161-input "
            "tower are not compatible. Retrain the model with the updated training "
            "code and use that new checkpoint for predict/perm."
        ) from exc
    model.eval()
    return model


def permute_feature_group(tensors, feature_name, rng):
    out = [x.clone() for x in tensors]
    feature_indices = [FEATURE_ORDER.index(name) for name in FEATURE_GROUPS[feature_name]]
    n = out[feature_indices[0]].shape[0]
    if n <= 1:
        return out

    perm = torch.from_numpy(rng.permutation(n)).long()
    for feature_idx in feature_indices:
        out[feature_idx] = out[feature_idx][perm]
    return out


def predict_one_chrom_with_permutation(
    case_dir,
    control_dir,
    variant_path,
    model,
    device,
    feature_name,
    permuted_sample,
    rng,
):
    variant_df = load_variant_table(variant_path)
    key_to_row = {k: i for i, k in enumerate(variant_df["key"].tolist())}

    control_shards = shard_files(control_dir)
    case_shards = shard_files(case_dir)
    if len(case_shards) != len(control_shards):
        raise ValueError(
            f"Different shard counts: case={len(case_shards)} control={len(control_shards)} "
            f"for {case_dir} vs {control_dir}"
        )

    result_parts = []
    for case_shard, control_shard in zip(case_shards, control_shards):
        case_keys, _, case_tensors = load_shard_records(case_shard, key_to_row=key_to_row)
        control_keys, _, control_tensors = load_shard_records(control_shard, key_to_row=key_to_row)

        if case_tensors is None:
            continue
        if control_tensors is None:
            raise ValueError(f"No matching control records for case shard {case_shard}")
        if list(case_keys) != list(control_keys):
            raise ValueError(
                "Case/control shard keys are not in the same order. "
                f"First case key={case_keys[0] if len(case_keys) else 'NA'}, "
                f"first control key={control_keys[0] if len(control_keys) else 'NA'}"
            )

        if permuted_sample in ("case", "both"):
            case_tensors = permute_feature_group(case_tensors, feature_name, rng)
        if permuted_sample in ("control", "both"):
            control_tensors = permute_feature_group(control_tensors, feature_name, rng)

        probs, pred_class = predict_feature_batch(case_tensors, control_tensors, model, device)
        row_idx = [key_to_row[k] for k in case_keys]
        rows = variant_df.iloc[row_idx, :].copy()
        result_parts.append(annotate_variant_rows(rows, probs, pred_class))

        del case_tensors, control_tensors
        gc.collect()

    if not result_parts:
        return pd.DataFrame()

    out = pd.concat(result_parts, axis=0, ignore_index=True)
    out["chrom"] = out["chrom"].astype(str)
    out["pos"] = out["pos"].astype(int)
    out = out.drop(columns=["key"], errors="ignore")
    return out


def run_permuted_predict(
    input_file_case,
    input_file_control,
    variant_file,
    model_path,
    output_file,
    feature_name,
    bam_path=None,
    target_chr=None,
    permuted_sample="both",
    device="cpu",
    seed=2026,
):
    check_feature_name(feature_name)
    target_chr = normalize_target_chr(target_chr)
    contigs = get_contigs_from_bam_or_default(bam_path)

    device = torch.device(device)
    rng = np.random.default_rng(seed)
    model = load_model(model_path, device)

    variant_df_list = []
    for chrom_id in target_chr:
        chrom = f"chr{chrom_id}"
        case_dir = input_file_case + f".{chrom}"
        control_dir = input_file_control + f".{chrom}"
        variant_path = variant_file + f".{chrom}"

        if not os.path.isdir(case_dir):
            print(f"Missing case shard directory {case_dir}; skip.")
            continue
        if not os.path.isdir(control_dir):
            print(f"Missing control shard directory {control_dir}; skip.")
            continue
        if not os.path.exists(variant_path):
            print(f"Missing variant file {variant_path}; skip.")
            continue

        print(f"Predicting {chrom} after permuting {feature_name} in {permuted_sample} sample(s)")
        pred_df = predict_one_chrom_with_permutation(
            case_dir=case_dir,
            control_dir=control_dir,
            variant_path=variant_path,
            model=model,
            device=device,
            feature_name=feature_name,
            permuted_sample=permuted_sample,
            rng=rng,
        )
        if pred_df.shape[0] == 0:
            print(f"No variants predicted for {chrom}; skip.")
            continue
        variant_df_list.append(pred_df)

    if not variant_df_list:
        raise RuntimeError("No variants were processed.")

    variant_df = pd.concat(variant_df_list, axis=0, ignore_index=True)
    variant_df.to_csv(output_file, index=False, sep="\t")

    filename = output_file[:-4] if output_file.endswith(".txt") else output_file
    variant_df_somatic = variant_df.loc[variant_df["Result_class"] == 3, :]
    variant_df_germline = variant_df.loc[variant_df["Result_class"].isin([1, 2]), :]

    variant_df_somatic.to_csv(filename + ".somatic.txt", index=False, sep="\t")
    variant_df_germline.to_csv(filename + ".germline.txt", index=False, sep="\t")

    convert_to_vcf(variant_df, filename + ".vcf", contigs)
    convert_to_vcf(variant_df_somatic, filename + ".somatic.vcf", contigs)
    convert_to_vcf(variant_df_germline, filename + ".germline.vcf", contigs)

    print("Permuted prediction details saved to", output_file)
    print("All permuted predictions VCF saved to", filename + ".vcf")
    print("Somatic permuted VCF saved to", filename + ".somatic.vcf")
    print("Germline permuted VCF saved to", filename + ".germline.vcf")
