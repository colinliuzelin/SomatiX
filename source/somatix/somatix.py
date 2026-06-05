#!/usr/bin/env python3

import argparse
import glob
import logging
import os
import shutil
import sys
from contextlib import redirect_stdout

try:
    from somatix.version import __version__
except Exception:
    __version__ = "unknown"


class SomatiXHelpFormatter(argparse.RawDescriptionHelpFormatter):
    """Show useful defaults while preserving multiline help text."""

    def _get_help_string(self, action):
        help_text = action.help
        if help_text is None:
            help_text = ""
        if "%(default)" not in help_text:
            has_default = action.default is not None and action.default is not argparse.SUPPRESS
            if action.option_strings and has_default:
                help_text += " (default: %(default)s)"
        return help_text


DEFAULT_THREADS = 8
DEFAULT_ALLELE_COUNTER = "allele_counter"
DEFAULT_MIN_ALT = 3
DEFAULT_MIN_VAF = 0.05
DEFAULT_MIN_TOTAL_CANDIDATE_COVERAGE = 3
DEFAULT_MIN_TOTAL_FEATURE_COVERAGE = 3
DEFAULT_MIN_MAPQ = 20
DEFAULT_MIN_BASEQ = 10
DEFAULT_MAX_DEPTH = 5000
DEFAULT_EXCL_FLAGS = 2316
DEFAULT_BLOCK_SIZE = 1_000_000
DEFAULT_FEATURE_DEPTH = 1000
DEFAULT_CHUNK_BP = 1000
DEFAULT_SHARD_SIZE = 100000
DEFAULT_COMPRESSION = None
DEFAULT_DEVICE = "cpu"
DEFAULT_PERM_SAMPLE = "both"
DEFAULT_PERM_SEED = 2026
WHOLE_GENOME_HELP = "Optional genomic interval, for example chr1:100000-200000. (default: whole genome)"


PERM_FEATURE_GROUPS = {
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


PERM_FEATURE_GROUP_DESCRIPTIONS = {
    "base_fraction": (
        "Five-channel base-fraction branch input: padded sequence context plus "
        "ref/alt A/C/G/T/DEL/INS base-fraction pileups on forward and reverse strands."
    ),
    "target_position_coverage": "Normalized total coverage at the candidate position.",
    "coverage": "Ref/alt allele coverage pileups on forward and reverse strands.",
    "mapq": "Mean mapping quality for ref/alt-supporting reads on forward and reverse strands.",
    "baseq": "Mean base quality for ref/alt-supporting reads on forward and reverse strands.",
    "mismatch_rate": "Mismatch-rate summaries for ref/alt-supporting reads on forward and reverse strands.",
}


def feature_group_help_text():
    lines = [
        "Experimental feature-testing command, not the standard SomatiX prediction workflow.",
        "Available --feature groups:",
    ]
    for name, members in PERM_FEATURE_GROUPS.items():
        lines.append(f"  {name}: {PERM_FEATURE_GROUP_DESCRIPTIONS[name]}")
        lines.append(f"    datasets: {', '.join(members)}")
    return "\n".join(lines)


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)


def ensure_directory_exists(directory):
    if directory and not os.path.exists(directory):
        logging.info(f"Creating directory: {directory}")
        os.makedirs(directory, exist_ok=True)


def check_required_files(files):
    for path in files:
        if path and not os.path.exists(path):
            logging.error(f"Required file does not exist: {path}")
            sys.exit(1)


def resolve_executable(path_or_name, label):
    if not path_or_name:
        logging.error(f"Missing executable for {label}")
        sys.exit(1)
    if os.path.sep in path_or_name:
        return path_or_name
    resolved = shutil.which(path_or_name)
    if resolved:
        return resolved
    return path_or_name


def normalize_chrom_name(chrom):
    chrom = str(chrom)
    if chrom.startswith("chr"):
        return chrom[3:]
    return chrom


def parse_target_chr_from_region(region=None):
    if region:
        chrom = region.split(":")[0]
        return [normalize_chrom_name(chrom)]
    return [str(i) for i in range(1, 23)]


def run_extract_candidates(args):
    """Extract candidate SNVs from one BAM using the allele_counter backend."""
    from somatix.bam2tpileup_selfAC import process_bam

    ensure_directory_exists(os.path.dirname(args.output) or ".")
    args.allele_counter = resolve_executable(args.allele_counter, "allele_counter")
    check_required_files([args.bam, args.ref, args.allele_counter])

    logging.info("Extracting candidate variants")
    logging.info(f"BAM: {args.bam}")
    logging.info(f"Region: {args.region if args.region else 'whole genome'}")
    logging.info(f"Output: {args.output}")
    logging.info(
        "Candidate filters: "
        f"min_alt={args.min_alt}, "
        f"min_vaf={args.min_vaf}, "
        f"min_total_coverage={args.min_total_coverage}, "
        f"min_mapq={args.min_mapq}, "
        f"min_baseq={args.min_baseq}, "
        f"max_depth={args.max_depth}"
    )

    with open(args.output, "w") as out:
        with redirect_stdout(out):
            process_bam(
                bam_path=args.bam,
                ref_path=args.ref,
                allele_counter_path=args.allele_counter,
                region=args.region,
                min_alt=args.min_alt,
                min_vaf=args.min_vaf,
                min_total_coverage=args.min_total_coverage,
                threads=args.threads,
                block_size=args.block_size,
                min_mapq=args.min_mapq,
                min_baseq=args.min_baseq,
                max_depth=args.max_depth,
                excl_flags=args.excl_flags,
            )

    logging.info("Candidate extraction finished")


def split_and_extract_features(
    candidate_file,
    bam,
    ref,
    output_prefix,
    filtered_prefix,
    threads,
    feature_depth,
    min_vaf,
    min_total_coverage,
    min_alt,
    shard_size,
    compression,
    chunk_bp,
    region=None,
):
    """
    Filter candidate variants by VAF/depth/ALT count, then extract HDF5 shard features.

    Output layout:
        {output_prefix}.chr1/shard_000000.h5
        {output_prefix}.chr1/manifest.tsv
        {filtered_prefix}.chr1

    The same filtered candidate table should be used for tumor and normal feature
    extraction so both HDF5 shard sets have identical loci and order.
    """
    from somatix.snp_feature_extraction_DNA_shard_light import (
        H5ShardWriter,
        extract_VCF_feature,
        filter_vcf,
    )

    ensure_directory_exists(os.path.dirname(output_prefix) or ".")
    ensure_directory_exists(os.path.dirname(filtered_prefix) or ".")

    target_chr = parse_target_chr_from_region(region)

    logging.info("Filtering candidate variants")
    logging.info(
        "Feature-input filters: "
        f"min_vaf={min_vaf}, "
        f"min_total_coverage={min_total_coverage}, "
        f"min_alt={min_alt}, "
        f"feature_depth={feature_depth}, "
        f"chunk_bp={chunk_bp}"
    )

    filtered_variant = filter_vcf(
        vcf_file=candidate_file,
        filtered_vcf_file=filtered_prefix,
        min_vaf=min_vaf,
        min_total_coverage=min_total_coverage,
        alt_min=min_alt,
    )

    processed_chr = []

    for chrom_id in target_chr:
        chrom = f"chr{normalize_chrom_name(chrom_id)}"
        each_filtered = filtered_variant.loc[
            filtered_variant["chrom"].astype(str) == chrom, :
        ].copy()

        if each_filtered.shape[0] == 0:
            logging.info(f"No candidate variants for {chrom}; skip")
            continue

        each_variant_file = filtered_prefix + f".{chrom}"
        each_output_dir = output_prefix + f".{chrom}"

        if os.path.isdir(each_output_dir):
            shutil.rmtree(each_output_dir)
        ensure_directory_exists(each_output_dir)

        each_filtered.to_csv(each_variant_file, sep="\t", index=False)

        logging.info(
            f"Extracting light shard features for {chrom}: "
            f"{each_filtered.shape[0]} candidates; "
            f"chunk_bp={chunk_bp}"
        )

        shard_writer = H5ShardWriter(
            output_dir=each_output_dir,
            shard_size=shard_size,
            compression=compression,
        )

        extract_VCF_feature(
            each_variant_file,
            bam,
            ref,
            threads,
            feature_depth,
            shard_writer,
            chunk_bp=chunk_bp,
        )

        shard_writer.close()
        processed_chr.append(normalize_chrom_name(chrom))

    logging.info(f"Feature extraction finished. Processed chromosomes: {processed_chr}")
    return processed_chr


def run_extract_features(args):
    check_required_files([args.candidates, args.bam, args.ref])

    split_and_extract_features(
        candidate_file=args.candidates,
        bam=args.bam,
        ref=args.ref,
        output_prefix=args.output_prefix,
        filtered_prefix=args.filtered_prefix,
        threads=args.threads,
        feature_depth=args.feature_depth,
        min_vaf=args.min_vaf,
        min_total_coverage=args.min_total_coverage,
        min_alt=args.min_alt,
        shard_size=args.shard_size,
        compression=args.compression,
        chunk_bp=args.chunk_bp,
        region=args.region,
    )


def run_predict(args):
    from somatix.predict_DNN_somatic import run_predict_dnn_somatic

    target_chr = parse_target_chr_from_region(args.region)
    check_required_files([args.model])

    if args.region:
        chrom = target_chr[0]
        check_required_files([args.variant_prefix + f".chr{chrom}"])
        for feature_dir in [
            args.case_features_prefix + f".chr{chrom}",
            args.control_features_prefix + f".chr{chrom}",
        ]:
            if not os.path.isdir(feature_dir):
                logging.error(f"Required feature shard directory does not exist: {feature_dir}")
                sys.exit(1)

    logging.info("Running SomatiX prediction from HDF5 feature shards")
    run_predict_dnn_somatic(
        input_file_case=args.case_features_prefix,
        input_file_control=args.control_features_prefix,
        variant_file=args.variant_prefix,
        model_path=args.model,
        output_file=args.output,
        bam_path=args.bam,
        target_chr=target_chr,
        device=args.device,
    )
    logging.info("Prediction finished")


def run_perm(args):
    from somatix.perm_importance import run_permuted_predict

    target_chr = parse_target_chr_from_region(args.region)
    check_required_files([args.model])

    if args.region:
        chrom = target_chr[0]
        check_required_files([args.variant_prefix + f".chr{chrom}"])
        for feature_dir in [
            args.case_features_prefix + f".chr{chrom}",
            args.control_features_prefix + f".chr{chrom}",
        ]:
            if not os.path.isdir(feature_dir):
                logging.error(f"Required feature shard directory does not exist: {feature_dir}")
                sys.exit(1)

    logging.info("Running SomatiX prediction with feature permutation")
    run_permuted_predict(
        input_file_case=args.case_features_prefix,
        input_file_control=args.control_features_prefix,
        variant_file=args.variant_prefix,
        model_path=args.model,
        output_file=args.output,
        feature_name=args.feature,
        bam_path=args.bam,
        target_chr=target_chr,
        permuted_sample=args.sample,
        device=args.device,
        seed=args.seed,
    )
    logging.info("Permuted prediction finished")


def run_call(args):
    ensure_directory_exists(args.outdir)
    args.allele_counter = resolve_executable(args.allele_counter, "allele_counter")

    required_files = [
        args.bam_case,
        args.bam_control,
        args.ref,
        args.allele_counter,
    ]

    if not args.skip_prediction:
        if not args.model:
            logging.error("--model is required unless --skip-prediction is used.")
            sys.exit(1)
        required_files.append(args.model)

    check_required_files(required_files)

    case_candidates = os.path.join(args.outdir, "case_candidates.txt")
    case_filtered_prefix = os.path.join(args.outdir, "case_candidates.filtered")
    case_features_prefix = os.path.join(args.outdir, "case_features")
    control_features_prefix = os.path.join(args.outdir, "control_features")
    predict_output = os.path.join(args.outdir, "somatix_predict.txt")

    logging.info("Step 1/4: extracting candidates from case BAM")
    candidate_args = argparse.Namespace(
        bam=args.bam_case,
        ref=args.ref,
        allele_counter=args.allele_counter,
        output=case_candidates,
        region=args.region,
        min_alt=args.min_alt,
        min_vaf=args.min_vaf,
        min_total_coverage=args.min_total_coverage,
        threads=args.threads,
        block_size=args.block_size,
        min_mapq=args.min_mapq,
        min_baseq=args.min_baseq,
        max_depth=args.max_depth,
        excl_flags=args.excl_flags,
    )
    run_extract_candidates(candidate_args)

    logging.info("Step 2/4: extracting case light shard features")
    processed_case_chr = split_and_extract_features(
        candidate_file=case_candidates,
        bam=args.bam_case,
        ref=args.ref,
        output_prefix=case_features_prefix,
        filtered_prefix=case_filtered_prefix,
        threads=args.threads,
        feature_depth=args.feature_depth,
        min_vaf=args.min_vaf,
        min_total_coverage=args.min_total_coverage,
        min_alt=args.min_alt,
        shard_size=args.shard_size,
        compression=args.compression,
        chunk_bp=args.chunk_bp,
        region=args.region,
    )

    logging.info("Step 3/4: extracting control light shard features at case candidate sites")
    # Important: the same case candidate table is used here so control features are
    # extracted at exactly the same candidate loci/order as the case features.
    processed_control_chr = split_and_extract_features(
        candidate_file=case_candidates,
        bam=args.bam_control,
        ref=args.ref,
        output_prefix=control_features_prefix,
        filtered_prefix=case_filtered_prefix,
        threads=args.threads,
        feature_depth=args.feature_depth,
        min_vaf=args.min_vaf,
        min_total_coverage=args.min_total_coverage,
        min_alt=args.min_alt,
        shard_size=args.shard_size,
        compression=args.compression,
        chunk_bp=args.chunk_bp,
        region=args.region,
    )

    final_chr = sorted(
        set(processed_case_chr).intersection(set(processed_control_chr)),
        key=lambda x: int(x) if str(x).isdigit() else 999,
    )
    if len(final_chr) == 0:
        logging.error("No chromosomes have both case and control features.")
        sys.exit(1)

    if args.skip_prediction:
        logging.info("Skipping Step 4/4: DL prediction was disabled by --skip-prediction")
        logging.info("SomatiX feature-generation finished")
        logging.info(f"Candidate file: {case_candidates}")
        logging.info(f"Filtered candidate prefix: {case_filtered_prefix}")
        logging.info(f"Case feature prefix: {case_features_prefix}")
        logging.info(f"Control feature prefix: {control_features_prefix}")
        logging.info(f"Processed chromosomes: {final_chr}")
        return

    logging.info(f"Step 4/4: predicting somatic variants for chromosomes: {final_chr}")
    from somatix.predict_DNN_somatic import run_predict_dnn_somatic

    run_predict_dnn_somatic(
        input_file_case=case_features_prefix,
        input_file_control=control_features_prefix,
        variant_file=case_filtered_prefix,
        model_path=args.model,
        output_file=predict_output,
        bam_path=args.bam_case,
        target_chr=final_chr,
        device=args.device,
    )

    logging.info("SomatiX call finished")
    logging.info(f"Final output prefix: {predict_output}")


def delete_related_files(outdir):
    patterns = [
        "case_candidates.txt",
        "case_candidates.filtered*",
        "case_features.chr*",
        "control_features.chr*",
        "*.tmp",
        "*.temp",
        "*.intermediate",
    ]

    for pattern in patterns:
        for path in glob.glob(os.path.join(outdir, pattern)):
            try:
                if os.path.isdir(path):
                    shutil.rmtree(path)
                else:
                    os.remove(path)
                logging.info(f"Deleted: {path}")
            except Exception as e:
                logging.error(f"Failed to delete {path}: {e}")


def add_common_candidate_args(parser):
    parser.add_argument(
        "--min-alt",
        type=int,
        default=DEFAULT_MIN_ALT,
        help="Minimum ALT-read count required for a candidate SNV.",
    )
    parser.add_argument(
        "--min-vaf",
        "--vaf",
        dest="min_vaf",
        type=float,
        default=DEFAULT_MIN_VAF,
        help="Minimum variant allele fraction required for a candidate SNV.",
    )
    parser.add_argument(
        "--min-total-coverage",
        "--min-total",
        dest="min_total_coverage",
        type=int,
        default=DEFAULT_MIN_TOTAL_CANDIDATE_COVERAGE,
        help="Minimum total read coverage required for a candidate SNV.",
    )
    parser.add_argument(
        "--min-mapq",
        "--min-MQ",
        type=int,
        default=DEFAULT_MIN_MAPQ,
        help="Minimum read mapping quality used during candidate extraction.",
    )
    parser.add_argument(
        "--min-baseq",
        "--min-BQ",
        type=int,
        default=DEFAULT_MIN_BASEQ,
        help="Minimum base quality used during candidate extraction.",
    )
    parser.add_argument(
        "--max-depth",
        type=int,
        default=DEFAULT_MAX_DEPTH,
        help="Maximum depth used during candidate extraction."
    )
    parser.add_argument(
        "--excl-flags",
        type=int,
        default=DEFAULT_EXCL_FLAGS,
        help="SAM flag bitmask for reads to exclude.",
    )
    parser.add_argument(
        "--block-size",
        type=int,
        default=DEFAULT_BLOCK_SIZE,
        help="Genomic block size used while scanning the BAM.",
    )


def add_common_feature_args(parser):
    parser.add_argument(
        "--feature-depth",
        "--depth",
        dest="feature_depth",
        type=int,
        default=DEFAULT_FEATURE_DEPTH,
        help="Maximum reads sampled per candidate site during feature extraction."
    )
    parser.add_argument(
        "--chunk-bp",
        type=int,
        default=DEFAULT_CHUNK_BP,
        help="Genomic chunk size for chunk-based feature extraction."
    )
    parser.add_argument(
        "--min-vaf",
        "--vaf",
        dest="min_vaf",
        type=float,
        default=DEFAULT_MIN_VAF,
        help="Minimum VAF retained when filtering the candidate table.",
    )
    parser.add_argument(
        "--min-total-coverage",
        "--min-total",
        dest="min_total_coverage",
        type=int,
        default=DEFAULT_MIN_TOTAL_FEATURE_COVERAGE,
        help="Minimum total coverage retained when filtering the candidate table.",
    )
    parser.add_argument(
        "--min-alt",
        type=int,
        default=DEFAULT_MIN_ALT,
        help="Minimum ALT-read count retained when filtering the candidate table.",
    )
    parser.add_argument(
        "--shard-size",
        type=int,
        default=DEFAULT_SHARD_SIZE,
        help="Maximum number of candidate records per HDF5 shard.",
    )
    parser.add_argument(
        "--compression",
        default=DEFAULT_COMPRESSION,
        choices=["lzf", "gzip"],
        help="HDF5 compression method for feature shards; omit to disable compression. (default: no compression)",
    )


def add_common_candidate_feature_args(parser):
    parser.add_argument(
        "--min-alt",
        type=int,
        default=DEFAULT_MIN_ALT,
        help="Minimum ALT-read count for candidate extraction and feature-input filtering.",
    )
    parser.add_argument(
        "--min-vaf",
        "--vaf",
        dest="min_vaf",
        type=float,
        default=DEFAULT_MIN_VAF,
        help="Minimum VAF for candidate extraction and feature-input filtering.",
    )
    parser.add_argument(
        "--min-total-coverage",
        "--min-total",
        dest="min_total_coverage",
        type=int,
        default=DEFAULT_MIN_TOTAL_CANDIDATE_COVERAGE,
        help="Minimum total coverage for candidate extraction and feature-input filtering.",
    )
    parser.add_argument(
        "--min-mapq",
        "--min-MQ",
        type=int,
        default=DEFAULT_MIN_MAPQ,
        help="Minimum read mapping quality used during candidate extraction.",
    )
    parser.add_argument(
        "--min-baseq",
        "--min-BQ",
        type=int,
        default=DEFAULT_MIN_BASEQ,
        help="Minimum base quality used during candidate extraction.",
    )
    parser.add_argument(
        "--max-depth",
        type=int,
        default=DEFAULT_MAX_DEPTH,
        help="Maximum depth used during candidate extraction."
    )
    parser.add_argument(
        "--excl-flags",
        type=int,
        default=DEFAULT_EXCL_FLAGS,
        help="SAM flag bitmask for reads to exclude.",
    )
    parser.add_argument(
        "--block-size",
        type=int,
        default=DEFAULT_BLOCK_SIZE,
        help="Genomic block size used while scanning the BAM.",
    )

    parser.add_argument(
        "--feature-depth",
        "--depth",
        dest="feature_depth",
        type=int,
        default=DEFAULT_FEATURE_DEPTH,
        help="Maximum reads sampled per candidate site during feature extraction."
    )
    parser.add_argument(
        "--chunk-bp",
        type=int,
        default=DEFAULT_CHUNK_BP,
        help="Genomic chunk size for chunk-based feature extraction."
    )

    parser.add_argument(
        "--shard-size",
        type=int,
        default=DEFAULT_SHARD_SIZE,
        help="Maximum number of candidate records per HDF5 shard.",
    )
    parser.add_argument(
        "--compression",
        default=DEFAULT_COMPRESSION,
        choices=["lzf", "gzip"],
        help="HDF5 compression method for feature shards; omit to disable compression. (default: no compression)",
    )


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog="somatix",
        description="SomatiX: platform-agnostic somatic SNV calling with tumor-control deep learning.",
        formatter_class=SomatiXHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"SomatiX {__version__}")

    subparsers = parser.add_subparsers(dest="command", required=True)

    p = subparsers.add_parser(
        "candidates",
        help="Extract candidate SNVs from BAM",
        formatter_class=SomatiXHelpFormatter,
    )
    p.add_argument("--bam", required=True, help="Input BAM file to scan for candidate SNVs.")
    p.add_argument("--ref", required=True, help="Indexed reference FASTA file.")
    p.add_argument("--output", required=True, help="Output tab-delimited candidate table.")
    p.add_argument("--region", help=WHOLE_GENOME_HELP)
    p.add_argument("--threads", type=int, default=DEFAULT_THREADS, help="Number of threads for candidate extraction.")
    p.add_argument("--allele-counter", default=DEFAULT_ALLELE_COUNTER, help="Path to the allele_counter executable.")
    add_common_candidate_args(p)

    p = subparsers.add_parser(
        "features",
        help="Extract reduced HDF5 shard features from candidate file",
        formatter_class=SomatiXHelpFormatter,
    )
    p.add_argument("--candidates", required=True, help="Input candidate table, usually from somatix candidates.")
    p.add_argument("--bam", required=True, help="BAM file used to extract features at candidate loci.")
    p.add_argument("--ref", required=True, help="Indexed reference FASTA file.")
    p.add_argument("--output-prefix", required=True, help="Prefix for per-chromosome HDF5 shard directories.")
    p.add_argument("--filtered-prefix", required=True, help="Prefix for filtered per-chromosome candidate tables.")
    p.add_argument("--threads", type=int, default=DEFAULT_THREADS, help="Number of threads for feature extraction.")
    p.add_argument("--region", help=WHOLE_GENOME_HELP)
    add_common_feature_args(p)

    p = subparsers.add_parser(
        "predict",
        help="Predict somatic variants from HDF5 shard features",
        formatter_class=SomatiXHelpFormatter,
    )
    p.add_argument("--case-features-prefix", required=True, help="Prefix for case/tumor feature shard directories.")
    p.add_argument("--control-features-prefix", required=True, help="Prefix for matched-normal/control feature shard directories.")
    p.add_argument("--variant-prefix", required=True, help="Prefix for filtered candidate tables used to align variants with feature shards.")
    p.add_argument("--model", required=True, help="PyTorch model checkpoint used for prediction.")
    p.add_argument("--output", required=True, help="Output prediction table path.")
    p.add_argument("--bam", help="Optional BAM path used to annotate VCF contigs.")
    p.add_argument("--device", default=DEFAULT_DEVICE, help="Device for model inference, such as cpu or cuda:0.")
    p.add_argument("--region", help="Optional genomic interval; prediction checks and runs only the selected chromosome. (default: whole genome)")

    p = subparsers.add_parser(
        "perm",
        help="Predict after permuting one selected HDF5 feature across samples",
        formatter_class=SomatiXHelpFormatter,
        epilog=feature_group_help_text(),
    )
    p.add_argument("--case-features-prefix", required=True, help="Prefix for case/tumor feature shard directories.")
    p.add_argument("--control-features-prefix", required=True, help="Prefix for matched-normal/control feature shard directories.")
    p.add_argument("--variant-prefix", required=True, help="Prefix for filtered candidate tables used to align variants with feature shards.")
    p.add_argument("--model", required=True, help="PyTorch model checkpoint used for prediction.")
    p.add_argument("--output", required=True, help="Output prediction table path.")
    p.add_argument(
        "--feature",
        required=True,
        choices=list(PERM_FEATURE_GROUPS.keys()),
        help="One feature group to permute across samples before prediction.",
    )
    p.add_argument(
        "--sample",
        choices=["case", "control", "both"],
        default=DEFAULT_PERM_SAMPLE,
        help="Which sample tower to permute for each selected feature.",
    )
    p.add_argument("--bam", help="Optional BAM path used to annotate VCF contigs.")
    p.add_argument("--device", default=DEFAULT_DEVICE, help="Device for model inference, such as cpu or cuda:0.")
    p.add_argument("--seed", type=int, default=DEFAULT_PERM_SEED, help="Random seed used to shuffle sample order.")
    p.add_argument("--region", help="Optional genomic interval; prediction checks and runs only the selected chromosome. (default: whole genome)")

    p = subparsers.add_parser(
        "call",
        help="Run full SomatiX pipeline",
        formatter_class=SomatiXHelpFormatter,
    )
    p.add_argument("--bam-case", required=True, help="Case/tumor BAM file.")
    p.add_argument("--bam-control", required=True, help="Matched normal/control BAM file.")
    p.add_argument("--ref", required=True, help="Indexed reference FASTA file.")
    p.add_argument("--model", required=False, help="PyTorch model checkpoint used for prediction. Required unless --skip-prediction is set.")
    p.add_argument("--outdir", required=True, help="Output directory for candidates, features and prediction results.")
    p.add_argument("--region", help=WHOLE_GENOME_HELP)
    p.add_argument("--threads", type=int, default=DEFAULT_THREADS, help="Number of threads used by candidate and feature extraction.")
    p.add_argument("--allele-counter", default=DEFAULT_ALLELE_COUNTER, help="Path to the allele_counter executable.")
    p.add_argument("--device", default=DEFAULT_DEVICE, help="Device for model inference, such as cpu or cuda:0.")
    p.add_argument(
        "--skip-prediction",
        action="store_true",
        help="Stop after candidate extraction and case/control feature extraction; do not run DL prediction.",
    )
    add_common_candidate_feature_args(p)

    p = subparsers.add_parser(
        "clean",
        help="Clean intermediate files",
        formatter_class=SomatiXHelpFormatter,
    )
    p.add_argument("--outdir", required=True, help="SomatiX output directory containing intermediate files to remove.")

    return parser.parse_args()


def main():
    args = parse_arguments()

    if args.command == "candidates":
        run_extract_candidates(args)
    elif args.command == "features":
        run_extract_features(args)
    elif args.command == "predict":
        run_predict(args)
    elif args.command == "perm":
        run_perm(args)
    elif args.command == "call":
        run_call(args)
    elif args.command == "clean":
        delete_related_files(args.outdir)
    else:
        raise ValueError(f"Unknown command: {args.command}")


if __name__ == "__main__":
    main()
