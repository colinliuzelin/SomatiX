import argparse
import gc
import glob
import math
import os
import warnings

import h5py
import numpy as np
import pandas as pd
import pysam

import torch
import torch.nn as nn
import torch.nn.functional as F


warnings.filterwarnings("ignore")


FEATURE_ORDER = [
    "one_hot_encoded_sequence",
    "target_position_coverage",

    "ref_base_fractions_forward",
    "ref_base_fractions_reverse",
    "alt_base_fractions_forward",
    "alt_base_fractions_reverse",

    "ref_coverage_forward",
    "ref_coverage_reverse",
    "alt_coverage_forward",
    "alt_coverage_reverse",

    "ref_mapq_forward",
    "ref_mapq_reverse",
    "alt_mapq_forward",
    "alt_mapq_reverse",

    "ref_baseq_forward",
    "ref_baseq_reverse",
    "alt_baseq_forward",
    "alt_baseq_reverse",

    "ref_mismatch_rate_forward",
    "ref_mismatch_rate_reverse",
    "alt_mismatch_rate_forward",
    "alt_mismatch_rate_reverse",


    "ref_target_coverage_forward",
    "ref_target_coverage_reverse",
    "alt_target_coverage_forward",
    "alt_target_coverage_reverse",
    
]


KEY_COLS = ["chrom", "pos", "ref", "alt"]


def make_key(chrom, pos, ref, alt):
    return f"{chrom}:{int(pos)}:{str(ref).upper()}:{str(alt).upper()}"


def decode_str_array(x):
    return np.array([v.decode("utf-8") if isinstance(v, bytes) else str(v) for v in x])


def normalize_target_chr(target_chr):
    if target_chr is None:
        return [str(i) for i in range(1, 23)]
    out = []
    for chrom in target_chr:
        s = str(chrom)
        if s.startswith("chr"):
            s = s[3:]
        out.append(s)
    return out


def get_contigs_from_bam_or_default(bam_path=None):
    if bam_path:
        bam = pysam.AlignmentFile(bam_path, "rb")
        contigs = {ref["SN"]: ref["LN"] for ref in bam.header["SQ"]}
        bam.close()
        return contigs

    return {
        "chr1": 248956422, "chr2": 242193529, "chr3": 198295559,
        "chr4": 190214555, "chr5": 181538259, "chr6": 170805979,
        "chr7": 159345973, "chr8": 145138636, "chr9": 138394717,
        "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
        "chr13": 114364328, "chr14": 107043718, "chr15": 101991189,
        "chr16": 90338345, "chr17": 83257441, "chr18": 80373285,
        "chr19": 58617616, "chr20": 64444167, "chr21": 46709983,
        "chr22": 50818468, "chrX": 156040895, "chrY": 57227415,
    }


def get_numeric_value(row, primary_name, fallback_names=None, default=0):
    """Return a numeric value from a pandas row using the current name first.

    This helper keeps VCF output robust to older intermediate files while
    preferring the current SomatiX candidate-table terminology.
    """
    if fallback_names is None:
        fallback_names = []

    for name in [primary_name] + list(fallback_names):
        if name in row.index and not pd.isna(row[name]):
            return row[name]

    return default


def convert_to_vcf(df, output_file, contigs):
    """
    Class definition:
        0 = 0/0 reference/artifact, FILTER=RefCall
        1 = 0/1 heterozygous germline, FILTER=GERMLINE
        2 = 1/1 homozygous germline, FILTER=GERMLINE
        3 = 1/1 somatic, FILTER=PASS

    VCF depth/frequency fields follow the updated candidate table:
        ref_reads, alt_reads, total_coverage, vaf
    """
    with open(output_file, "w") as vcf:
        vcf.write("##fileformat=VCFv4.2\n")
        vcf.write("##source=SomatiX\n")

        for chrom, length in contigs.items():
            vcf.write(f"##contig=<ID={chrom},length={length}>\n")

        vcf.write("##FILTER=<ID=PASS,Description=\"SomatiX somatic call.\">\n")
        vcf.write("##FILTER=<ID=RefCall,Description=\"SomatiX classifies this candidate as reference/artifact.\">\n")
        vcf.write("##FILTER=<ID=NoCall,Description=\"Candidate has zero total coverage.\">\n")
        vcf.write("##FILTER=<ID=GERMLINE,Description=\"SomatiX classifies this candidate as germline.\">\n")

        vcf.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total candidate-site read depth from allele counting.\">\n")
        vcf.write("##INFO=<ID=VAF,Number=A,Type=Float,Description=\"Variant allele fraction, ALT reads divided by total coverage.\">\n")
        vcf.write("##INFO=<ID=REF_READS,Number=1,Type=Integer,Description=\"Number of reads supporting the reference allele.\">\n")
        vcf.write("##INFO=<ID=ALT_READS,Number=1,Type=Integer,Description=\"Number of reads supporting the selected ALT allele.\">\n")
        vcf.write("##INFO=<ID=C0,Number=1,Type=Float,Description=\"SomatiX probability of class 0: reference/artifact.\">\n")
        vcf.write("##INFO=<ID=C1,Number=1,Type=Float,Description=\"SomatiX probability of class 1: heterozygous germline.\">\n")
        vcf.write("##INFO=<ID=C2,Number=1,Type=Float,Description=\"SomatiX probability of class 2: homozygous germline.\">\n")
        vcf.write("##INFO=<ID=C3,Number=1,Type=Float,Description=\"SomatiX probability of class 3: somatic.\">\n")

        vcf.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"SomatiX genotype-like class encoding.\">\n")
        vcf.write("##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"Phred-scaled confidence of the predicted SomatiX class.\">\n")
        vcf.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total candidate-site read depth from allele counting.\">\n")
        vcf.write("##FORMAT=<ID=AD,Number=2,Type=Integer,Description=\"Allelic depths for REF and selected ALT.\">\n")
        vcf.write("##FORMAT=<ID=VAF,Number=A,Type=Float,Description=\"Variant allele fraction, ALT reads divided by total coverage.\">\n")
        vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")

        for _, row in df.iterrows():
            chrom = row["chrom"]
            pos = int(row["pos"])
            ref = str(row["ref"]).upper()
            alt = str(row["alt"]).upper()
            result_class = int(row["Result_class"])
            probs = [float(row[f"class_{i}"]) for i in range(4)]

            max_prob = max(probs)
            error_prob = max(1e-10, 1.0 - max_prob)
            qual = min(-10.0 * math.log10(error_prob), 99.0)
            gq = qual

            ref_reads = int(get_numeric_value(row, "ref_reads", default=0))
            alt_reads = int(get_numeric_value(row, "alt_reads", default=0))
            dp = int(get_numeric_value(
                row,
                "total_coverage",
                fallback_names=["DP", "depth"],
                default=ref_reads + alt_reads,
            ))

            vaf = float(get_numeric_value(
                row,
                "vaf",
                fallback_names=["VAF", "AF", "ratio_1"],
                default=(alt_reads / dp if dp > 0 else 0.0),
            ))

            ad = f"{ref_reads},{alt_reads}"

            if dp == 0:
                gt, filter_status = "./.", "NoCall"
            elif result_class == 0:
                gt, filter_status = "0/0", "RefCall"
            elif result_class == 1:
                gt, filter_status = "0/1", "GERMLINE"
            elif result_class == 2:
                gt, filter_status = "1/1", "GERMLINE"
            elif result_class == 3:
                gt, filter_status = "1/1", "PASS"
            else:
                gt, filter_status = "./.", "NoCall"

            info = (
                f"DP={dp};VAF={vaf:.4f};"
                f"REF_READS={ref_reads};ALT_READS={alt_reads};"
                f"C0={probs[0]:.4f};C1={probs[1]:.4f};"
                f"C2={probs[2]:.4f};C3={probs[3]:.4f}"
            )
            sample = f"{gt}:{gq:.2f}:{dp}:{ad}:{vaf:.4f}"
            vcf.write(
                f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t{qual:.2f}\t"
                f"{filter_status}\t{info}\tGT:GQ:DP:AD:VAF\t{sample}\n"
            )


# -----------------------------
# Model
# -----------------------------
class CombineTower(nn.Module):
    def __init__(self):
        super().__init__()

        self.allele_conv1 = nn.Conv2d(5, 64, kernel_size=(6, 3), padding=(0, 1))
        self.allele_pool1 = nn.MaxPool2d((1, 2))
        self.allele_conv2 = nn.Conv2d(64, 64, kernel_size=(1, 3), padding=(0, 1))
        self.allele_pool2 = nn.MaxPool2d((1, 2))

        self.dropout = nn.Dropout(0.3)

        self.ref_cov_fc1_fwd = nn.Linear(61, 32)
        self.ref_cov_fc1_rev = nn.Linear(61, 32)
        self.ref_cov_fc2 = nn.Linear(64, 16)

        self.alt_cov_fc1_fwd = nn.Linear(61, 32)
        self.alt_cov_fc1_rev = nn.Linear(61, 32)
        self.alt_cov_fc2 = nn.Linear(64, 16)

        self.mapq_fc1 = nn.Linear(4, 4)

        self.ref_baseq_fc1_fwd = nn.Linear(61, 32)
        self.ref_baseq_fc1_rev = nn.Linear(61, 32)
        self.ref_baseq_fc2 = nn.Linear(64, 16)

        self.alt_baseq_fc1_fwd = nn.Linear(61, 32)
        self.alt_baseq_fc1_rev = nn.Linear(61, 32)
        self.alt_baseq_fc2 = nn.Linear(64, 16)

        self.mismatch_fc1 = nn.Linear(4, 4)

        # total = 960 + 1 + 68 + 4 = 1033
        self.fc1 = nn.Linear(1033, 128)

    def forward(
        self,
        one_hot_seq,
        target_position_coverage,
        ref_frac_fwd,
        ref_frac_rev,
        alt_frac_fwd,
        alt_frac_rev,
        ref_cov_fwd,
        ref_cov_rev,
        alt_cov_fwd,
        alt_cov_rev,
        ref_mapq_fwd,
        ref_mapq_rev,
        alt_mapq_fwd,
        alt_mapq_rev,
        ref_baseq_fwd,
        ref_baseq_rev,
        alt_baseq_fwd,
        alt_baseq_rev,

        ref_mismatch_fwd,
        ref_mismatch_rev,
        alt_mismatch_fwd,
        alt_mismatch_rev,
    ):

        zero_pad = torch.zeros(
        (one_hot_seq.size(0), 2, one_hot_seq.size(2)),
        dtype=one_hot_seq.dtype,
        device=one_hot_seq.device,
        )
        
        one_hot_seq_6 = torch.cat([one_hot_seq, zero_pad], dim=1)
        
        allele_input = torch.stack(
            [
                one_hot_seq_6,
                ref_frac_fwd,
                ref_frac_rev,
                alt_frac_fwd,
                alt_frac_rev,
            ],
            dim=1,
        )
    
        x = F.relu(self.allele_conv1(allele_input))
        x = self.allele_pool1(x)
        x = self.dropout(x)

        x = F.relu(self.allele_conv2(x))
        x = self.allele_pool2(x)
        x = self.dropout(x)

        x = x.view(x.size(0), -1)

        target_position_coverage = target_position_coverage.view(
            target_position_coverage.size(0), -1
        )

        ref_cov = F.relu(self.ref_cov_fc2(torch.cat([
            F.relu(self.ref_cov_fc1_fwd(ref_cov_fwd)),
            F.relu(self.ref_cov_fc1_rev(ref_cov_rev)),
        ], dim=1)))

        alt_cov = F.relu(self.alt_cov_fc2(torch.cat([
            F.relu(self.alt_cov_fc1_fwd(alt_cov_fwd)),
            F.relu(self.alt_cov_fc1_rev(alt_cov_rev)),
        ], dim=1)))

        mapq = torch.cat([
            ref_mapq_fwd.view(ref_mapq_fwd.size(0), -1),
            ref_mapq_rev.view(ref_mapq_rev.size(0), -1),
            alt_mapq_fwd.view(alt_mapq_fwd.size(0), -1),
            alt_mapq_rev.view(alt_mapq_rev.size(0), -1),
        ], dim=1)
        mapq = F.relu(self.mapq_fc1(mapq))

        ref_baseq = F.relu(self.ref_baseq_fc2(torch.cat([
            F.relu(self.ref_baseq_fc1_fwd(ref_baseq_fwd)),
            F.relu(self.ref_baseq_fc1_rev(ref_baseq_rev)),
        ], dim=1)))

        alt_baseq = F.relu(self.alt_baseq_fc2(torch.cat([
            F.relu(self.alt_baseq_fc1_fwd(alt_baseq_fwd)),
            F.relu(self.alt_baseq_fc1_rev(alt_baseq_rev)),
        ], dim=1)))

        mismatch = torch.cat([
            ref_mismatch_fwd.view(ref_mismatch_fwd.size(0), -1),
            ref_mismatch_rev.view(ref_mismatch_rev.size(0), -1),
            alt_mismatch_fwd.view(alt_mismatch_fwd.size(0), -1),
            alt_mismatch_rev.view(alt_mismatch_rev.size(0), -1),
        ], dim=1)
        mismatch = F.relu(self.mismatch_fc1(mismatch))

        combined = torch.cat(
            (
                x,
                target_position_coverage,
                ref_cov,
                alt_cov,
                mapq,
                ref_baseq,
                alt_baseq,
                mismatch,
            ),
            dim=1,
        )

        return F.relu(self.fc1(combined))


class DNASeqCoverageModel_somatic(nn.Module):
    def __init__(self):
        super().__init__()

        self.tower_R = CombineTower()
        self.tower_N = CombineTower()

        # New tumor + normal target-position raw coverage branch:
        # 8 raw values -> normalize by per-example max -> 8 hidden -> dropout -> 4 output
        self.target_cov_fc1 = nn.Linear(8, 8)
        self.target_cov_fc2 = nn.Linear(8, 4)
        self.target_cov_dropout = nn.Dropout(0.05)

        # R_feat + N_feat + target_cov = 128 + 128= 256
        self.fc2 = nn.Linear(128 * 2 , 128)
        self.fc3 = nn.Linear(128, 64)
        self.fc_output = nn.Linear(64, 4)

    def forward(
        self,
        R_one_hot_seq,
        R_target_position_coverage,
        R_ref_frac_fwd,
        R_ref_frac_rev,
        R_alt_frac_fwd,
        R_alt_frac_rev,
        R_ref_cov_fwd,
        R_ref_cov_rev,
        R_alt_cov_fwd,
        R_alt_cov_rev,
        R_ref_mapq_fwd,
        R_ref_mapq_rev,
        R_alt_mapq_fwd,
        R_alt_mapq_rev,
        R_ref_baseq_fwd,
        R_ref_baseq_rev,
        R_alt_baseq_fwd,
        R_alt_baseq_rev,

        R_ref_mismatch_fwd,
        R_ref_mismatch_rev,
        R_alt_mismatch_fwd,
        R_alt_mismatch_rev,

        R_ref_target_cov_fwd,
        R_ref_target_cov_rev,
        R_alt_target_cov_fwd,
        R_alt_target_cov_rev,

        N_one_hot_seq,
        N_target_position_coverage,
        N_ref_frac_fwd,
        N_ref_frac_rev,
        N_alt_frac_fwd,
        N_alt_frac_rev,
        N_ref_cov_fwd,
        N_ref_cov_rev,
        N_alt_cov_fwd,
        N_alt_cov_rev,
        N_ref_mapq_fwd,
        N_ref_mapq_rev,
        N_alt_mapq_fwd,
        N_alt_mapq_rev,
        N_ref_baseq_fwd,
        N_ref_baseq_rev,
        N_alt_baseq_fwd,
        N_alt_baseq_rev,

        N_ref_mismatch_fwd,
        N_ref_mismatch_rev,
        N_alt_mismatch_fwd,
        N_alt_mismatch_rev,

        N_ref_target_cov_fwd,
        N_ref_target_cov_rev,
        N_alt_target_cov_fwd,
        N_alt_target_cov_rev,
    ):
        R_feat = self.tower_R(
            R_one_hot_seq,
            R_target_position_coverage,
            R_ref_frac_fwd,
            R_ref_frac_rev,
            R_alt_frac_fwd,
            R_alt_frac_rev,
            R_ref_cov_fwd,
            R_ref_cov_rev,
            R_alt_cov_fwd,
            R_alt_cov_rev,
            R_ref_mapq_fwd,
            R_ref_mapq_rev,
            R_alt_mapq_fwd,
            R_alt_mapq_rev,
            R_ref_baseq_fwd,
            R_ref_baseq_rev,
            R_alt_baseq_fwd,
            R_alt_baseq_rev,

            R_ref_mismatch_fwd,
            R_ref_mismatch_rev,
            R_alt_mismatch_fwd,
            R_alt_mismatch_rev,
        )

        N_feat = self.tower_N(
            N_one_hot_seq,
            N_target_position_coverage,
            N_ref_frac_fwd,
            N_ref_frac_rev,
            N_alt_frac_fwd,
            N_alt_frac_rev,
            N_ref_cov_fwd,
            N_ref_cov_rev,
            N_alt_cov_fwd,
            N_alt_cov_rev,
            N_ref_mapq_fwd,
            N_ref_mapq_rev,
            N_alt_mapq_fwd,
            N_alt_mapq_rev,
            N_ref_baseq_fwd,
            N_ref_baseq_rev,
            N_alt_baseq_fwd,
            N_alt_baseq_rev,

            N_ref_mismatch_fwd,
            N_ref_mismatch_rev,
            N_alt_mismatch_fwd,
            N_alt_mismatch_rev,
        )

        target_cov = torch.cat([
            R_ref_target_cov_fwd.view(R_ref_target_cov_fwd.size(0), -1),
            R_ref_target_cov_rev.view(R_ref_target_cov_rev.size(0), -1),
            R_alt_target_cov_fwd.view(R_alt_target_cov_fwd.size(0), -1),
            R_alt_target_cov_rev.view(R_alt_target_cov_rev.size(0), -1),
            N_ref_target_cov_fwd.view(N_ref_target_cov_fwd.size(0), -1),
            N_ref_target_cov_rev.view(N_ref_target_cov_rev.size(0), -1),
            N_alt_target_cov_fwd.view(N_alt_target_cov_fwd.size(0), -1),
            N_alt_target_cov_rev.view(N_alt_target_cov_rev.size(0), -1),
        ], dim=1)

        #max_target_cov = torch.max(target_cov, dim=1, keepdim=True)[0]
        #target_cov = target_cov / (max_target_cov + 1e-10)

        #target_cov = F.relu(self.target_cov_fc1(target_cov))
        #target_cov = self.target_cov_dropout(target_cov)
        #target_cov = F.relu(self.target_cov_fc2(target_cov))

        #x = torch.cat((R_feat, N_feat, target_cov), dim=1)

        x = torch.cat((R_feat, N_feat), dim=1)
        
        x = F.relu(self.fc2(x))
        x = F.dropout(x, p=0.1, training=self.training)
        
        x = F.relu(self.fc3(x))
        x = F.dropout(x, p=0.1, training=self.training)

        return self.fc_output(x)

def h5_to_tensor_batch(h5, indices=None):
    arrays = []
    for name in FEATURE_ORDER:
        data = h5[name][:] if indices is None else h5[name][indices]
        arrays.append(torch.from_numpy(data.astype(np.float32, copy=False)))
    return arrays


def load_shard_records(shard_path, key_to_row=None):
    with h5py.File(shard_path, "r") as h5:
        keys = decode_str_array(h5["key"][:])
        if key_to_row is None:
            keep = np.arange(len(keys))
        else:
            keep = np.array([i for i, k in enumerate(keys) if k in key_to_row], dtype=np.int64)
        if keep.size == 0:
            return keys, keep, None
        tensors = h5_to_tensor_batch(h5, keep)
        kept_keys = keys[keep]
    return kept_keys, keep, tensors


def read_manifest_keys(feature_dir):
    manifest = os.path.join(feature_dir, "manifest.tsv")
    if not os.path.exists(manifest):
        return None
    df = pd.read_csv(manifest, sep="\t")
    return df["key"].astype(str).tolist()


def shard_files(feature_dir):
    files = sorted(glob.glob(os.path.join(feature_dir, "shard_*.h5")))
    if not files:
        raise FileNotFoundError(f"No shard_*.h5 files found in {feature_dir}")
    return files


def load_variant_table(variant_path):
    df = pd.read_csv(variant_path, sep="\t")
    missing = [c for c in KEY_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"Variant file {variant_path} missing columns: {missing}")
    df = df.copy()
    df["key"] = [make_key(r.chrom, r.pos, r.ref, r.alt) for r in df.itertuples(index=False)]
    df = df.drop_duplicates(subset="key", keep="first").reset_index(drop=True)
    return df


def predict_feature_batch(case_tensors, control_tensors, model, device):
    inputs = [x.to(device) for x in case_tensors + control_tensors]
    with torch.no_grad():
        logits = model(*inputs)
        probs = torch.softmax(logits, dim=1).cpu().numpy()
    pred_class = np.argmax(probs, axis=1)
    return probs, pred_class


def annotate_variant_rows(rows, probs, pred_class):
    out = rows.copy().reset_index(drop=True)
    for i in range(4):
        out[f"class_{i}"] = probs[:, i]
    out["Result_reference"] = (pred_class == 0).astype(int)
    out["Result_hete_germline"] = (pred_class == 1).astype(int)
    out["Result_homo_germline"] = (pred_class == 2).astype(int)
    out["Result_germline"] = ((pred_class == 1) | (pred_class == 2)).astype(int)
    out["Result_somatic"] = (pred_class == 3).astype(int)
    out["Result_class"] = pred_class
    return out


def predict_one_chrom(case_dir, control_dir, variant_path, model, device):
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


def run_predict_dnn_somatic(
    input_file_case,
    input_file_control,
    variant_file,
    model_path,
    output_file,
    bam_path=None,
    target_chr=None,
    device="cpu",
):
    target_chr = normalize_target_chr(target_chr)
    contigs = get_contigs_from_bam_or_default(bam_path)

    device = torch.device(device)
    model = DNASeqCoverageModel_somatic().to(device)
    state = torch.load(model_path, map_location=device)
    if isinstance(state, dict) and "state_dict" in state:
        state = state["state_dict"]
    model.load_state_dict(state)
    model.eval()
    print("SomatiX total-fraction + MAPQ + baseQ model loaded for predictions.")

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

        print(f"Predicting {chrom}")
        pred_df = predict_one_chrom(case_dir, control_dir, variant_path, model, device)
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

    print("Prediction details saved to", output_file)
    print("All predictions VCF saved to", filename + ".vcf")
    print("Somatic VCF saved to", filename + ".somatic.vcf")
    print("Germline VCF saved to", filename + ".germline.vcf")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run SomatiX prediction from HDF5 feature shards.")
    parser.add_argument("--input_file_case", required=True, help="Case feature directory prefix, e.g. case_features")
    parser.add_argument("--input_file_control", required=True, help="Control feature directory prefix, e.g. control_features")
    parser.add_argument("--variant_file", required=True, help="Filtered variant file prefix")
    parser.add_argument("--model_path", required=True)
    parser.add_argument("--output_file", required=True)
    parser.add_argument("--bam", required=False)
    parser.add_argument("--device", default="cpu", help="cpu or cuda")
    parser.add_argument("--target_chr", nargs="+", default=[str(i) for i in range(1, 23)])
    args = parser.parse_args()

    run_predict_dnn_somatic(
        input_file_case=args.input_file_case,
        input_file_control=args.input_file_control,
        variant_file=args.variant_file,
        model_path=args.model_path,
        output_file=args.output_file,
        bam_path=args.bam,
        target_chr=args.target_chr,
        device=args.device,
    )
