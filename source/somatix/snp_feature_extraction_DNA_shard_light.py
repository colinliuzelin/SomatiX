import argparse
import os
import h5py
import pysam
import numpy as np
from multiprocessing import Pool
import pandas as pd
import progressbar
from bisect import bisect_left


GLOBAL_BAM = None
GLOBAL_REF = None


BASE_TO_IDX = {
    "A": 0,
    "C": 1,
    "G": 2,
    "T": 3,
    "DEL": 4,
    "INS": 5,
}


FEATURE_ORDER = [
    "one_hot_encoded_sequence",
    "total_base_fraction",
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

    "ref_read_length_forward",
    "ref_read_length_reverse",
    "alt_read_length_forward",
    "alt_read_length_reverse",

    "ref_target_relative_position_forward",
    "ref_target_relative_position_reverse",
    "alt_target_relative_position_forward",
    "alt_target_relative_position_reverse",

    "ref_target_coverage_forward",
    "ref_target_coverage_reverse",
    "alt_target_coverage_forward",
    "alt_target_coverage_reverse",
]


def init_worker(bam_file, reference_file):
    global GLOBAL_BAM, GLOBAL_REF
    GLOBAL_BAM = pysam.AlignmentFile(bam_file, "rb")
    GLOBAL_REF = pysam.FastaFile(reference_file)


def read_passes_basic_filter(read, min_mapq=20):
    if (
        read.is_unmapped or
        read.is_duplicate or
        read.mapping_quality < min_mapq or
        read.is_supplementary
    ):
        return False

    if read.query_sequence is None:
        return False

    if read.cigartuples is None:
        return False

    if read.reference_start is None or read.reference_end is None:
        return False

    return True


def filter_vcf(
    vcf_file,
    filtered_vcf_file,
    min_vaf,
    min_total_coverage,
    alt_min,
):
    vcf_data = pd.read_csv(vcf_file, sep="\t")

    filtered_vcf = vcf_data[
        (vcf_data["vaf"] >= min_vaf) &
        (vcf_data["alt_reads"] >= alt_min) &
        (vcf_data["total_coverage"] >= min_total_coverage) &
        (vcf_data["ref"].astype(str).str.len() == 1) &
        (vcf_data["alt"].astype(str).str.len() == 1)
    ]

    df_filtered = filtered_vcf[
        filtered_vcf["chrom"].astype(str).str.match(r"^chr([1-9]|1[0-9]|2[0-2])$")
    ].copy()

    chrom_order = {f"chr{i}": i for i in range(1, 23)}
    df_filtered["chrom_order"] = df_filtered["chrom"].map(chrom_order)

    df_sorted = (
        df_filtered
        .sort_values(by=["chrom_order", "pos"])
        .drop(columns="chrom_order")
        .reset_index(drop=True)
    )

    return df_sorted


def extract_variants_preserve_order(vcf_file):
    variants = []

    with open(vcf_file, "r") as f_vcf:
        _ = f_vcf.readline()
        line_index = 0

        for line in f_vcf:
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")
            if len(fields) < 4:
                continue

            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[2]
            alt = fields[3]

            variants.append({
                "index": line_index,
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
            })

            line_index += 1

    return variants


def make_chunk_tasks(vcf_file, bam_file, reference_file, read_depth, chunk_bp):
    """
    Build chunk-level tasks.

    Variants are grouped by:
      - chromosome
      - fixed genomic bin: (pos - 1) // chunk_bp

    Each chunk task fetches BAM reads once over the min/max candidate position
    range in that chunk, then reuses those reads for all variants in the chunk.
    """
    variants = extract_variants_preserve_order(vcf_file)

    if len(variants) == 0:
        return [], 0

    grouped = {}

    for v in variants:
        chrom = v["chrom"]
        bin_id = (int(v["pos"]) - 1) // int(chunk_bp)
        key = (chrom, bin_id)

        if key not in grouped:
            grouped[key] = []

        grouped[key].append(v)

    chunk_tasks = []

    for chunk_id, key in enumerate(sorted(grouped.keys(), key=lambda x: (x[0], x[1]))):
        chrom, bin_id = key
        chunk_variants = sorted(grouped[key], key=lambda x: x["pos"])

        min_pos = min(v["pos"] for v in chunk_variants)
        max_pos = max(v["pos"] for v in chunk_variants)

        # pysam.fetch uses 0-based half-open coordinates.
        # Fetching min_pos-1 to max_pos is sufficient to collect reads
        # overlapping candidate target positions.
        fetch_start = max(0, min_pos - 1)
        fetch_end = max_pos

        chunk_tasks.append((
            chunk_id,
            bam_file,
            reference_file,
            chrom,
            fetch_start,
            fetch_end,
            chunk_variants,
            read_depth,
        ))

    return chunk_tasks, len(variants)


def get_reference_sequence(chrom, pos, window=30):
    global GLOBAL_REF

    start = max(0, pos - 1 - window)
    end = pos + window
    return GLOBAL_REF.fetch(chrom, start, end).upper()


def one_hot_encode_sequence(sequence):
    mapping = {
        "A": [1, 0, 0, 0],
        "C": [0, 1, 0, 0],
        "G": [0, 0, 1, 0],
        "T": [0, 0, 0, 1],
    }

    one_hot_matrix = np.array(
        [mapping.get(base, [0, 0, 0, 0]) for base in sequence],
        dtype=np.float32
    )

    return one_hot_matrix.T


def init_feature_accumulators():
    return {
        "base_counts": np.zeros((6, 61), dtype=np.float32),
        "coverage": np.zeros(61, dtype=np.float32),
        "fraction_denominator": np.zeros(61, dtype=np.float32),
        "baseq_sum": np.zeros(61, dtype=np.float32),
        "baseq_count": np.zeros(61, dtype=np.float32),
        "mapq_sum": 0.0,
        "read_count": 0,
        "mismatch_rate_sum": 0.0,
        "read_length_sum": 0.0,
        "target_relative_position_sum": 0.0,
        "target_coverage_count": 0.0,
    }


def get_read_mismatch_rate(read):
    aligned_len = float(read.query_alignment_length or 0)
    if aligned_len <= 0:
        return 0.0

    try:
        nm = float(read.get_tag("NM"))
    except KeyError:
        return 0.0

    return min(nm / aligned_len, 1.0)


def get_read_length(read):
    read_len = read.query_length
    if read_len is None or read_len <= 0:
        return 0.0
    return float(read_len)


def get_target_relative_position(read, pos):
    read_len = get_read_length(read)
    if read_len <= 0 or read.reference_start is None:
        return 0.0

    rel_pos = float((pos - 1) - read.reference_start) / read_len
    return float(np.clip(rel_pos, 0.0, 1.0))


def add_read_to_accumulator(acc, read, window_items, pos):
    read_bases = read.query_sequence
    read_qualities = read.query_qualities

    acc["mapq_sum"] += float(read.mapping_quality)
    acc["read_count"] += 1
    acc["mismatch_rate_sum"] += get_read_mismatch_rate(read)
    acc["read_length_sum"] += get_read_length(read)
    acc["target_relative_position_sum"] += get_target_relative_position(read, pos)

    pos0 = pos - 1

    base_counts = acc["base_counts"]
    coverage = acc["coverage"]
    fraction_denominator = acc["fraction_denominator"]
    baseq_sum = acc["baseq_sum"]
    baseq_count = acc["baseq_count"]

    for window_idx, query_pos, ref_pos in window_items:
        if ref_pos == pos0 and query_pos is not None:
            acc["target_coverage_count"] += 1.0

        if query_pos is None:
            base_counts[BASE_TO_IDX["DEL"], window_idx] += 1.0
            fraction_denominator[window_idx] += 1.0
            continue

        if query_pos == "INS":
            base_counts[BASE_TO_IDX["INS"], window_idx] += 1.0
            fraction_denominator[window_idx] += 1.0
            continue

        if query_pos >= len(read_bases):
            continue

        base = read_bases[query_pos].upper()
        base_idx = BASE_TO_IDX.get(base)

        if base_idx is None:
            continue

        base_counts[base_idx, window_idx] += 1.0
        coverage[window_idx] += 1.0
        fraction_denominator[window_idx] += 1.0

        if read_qualities is not None and query_pos < len(read_qualities):
            baseq_sum[window_idx] += float(read_qualities[query_pos])
            baseq_count[window_idx] += 1.0


def finalize_base_fraction_and_coverage(acc):
    base_fractions = acc["base_counts"] / (
        acc["fraction_denominator"].reshape(1, -1) + 1e-10
    )
    return base_fractions.astype(np.float32), acc["coverage"].astype(np.float32)


def finalize_average_baseq(acc):
    mean_baseq = acc["baseq_sum"] / (acc["baseq_count"] + 1e-10)
    mean_baseq = np.minimum(mean_baseq / 40.0, 1.0)
    return mean_baseq.astype(np.float32)


def finalize_average_mapq(acc):
    if acc["read_count"] == 0:
        return np.array([0.0], dtype=np.float32)

    return np.array(
        [(acc["mapq_sum"] / float(acc["read_count"])) / 60.0],
        dtype=np.float32
    )


def finalize_average_mismatch_rate(acc):
    if acc["read_count"] == 0:
        return np.array([0.0], dtype=np.float32)

    return np.array(
        [acc["mismatch_rate_sum"] / float(acc["read_count"])],
        dtype=np.float32
    )


def finalize_average_read_length(acc, cap=20000.0):
    if acc["read_count"] == 0:
        return np.array([0.0], dtype=np.float32)

    mean_len = acc["read_length_sum"] / float(acc["read_count"])
    norm_len = min(mean_len, cap) / cap

    return np.array([norm_len], dtype=np.float32)


def finalize_average_target_relative_position(acc):
    if acc["read_count"] == 0:
        return np.array([0.0], dtype=np.float32)

    return np.array(
        [acc["target_relative_position_sum"] / float(acc["read_count"])],
        dtype=np.float32
    )


def finalize_target_coverage_count(acc):
    return np.array([acc["target_coverage_count"]], dtype=np.float32)


def extract_read_target_and_window(read, pos, window=30):
    """
    Fast CIGAR-based replacement for read.get_aligned_pairs().

    Feature behavior:
      - A/C/G/T are counted per aligned base.
      - DEL is counted per deleted/skipped reference base.
      - INS is counted per inserted base, not per insertion event.
      - INS is anchored to the previous reference position.
      - base_at_target is "INS" if insertion is anchored at the target.
      - base_at_target is "DEL" if target reference position has query_pos=None.
    """
    read_bases = read.query_sequence
    if read_bases is None:
        return None, None

    cigartuples = read.cigartuples
    if cigartuples is None:
        return None, None

    ref_start = read.reference_start
    if ref_start is None:
        return None, None

    pos0 = pos - 1
    window_start = pos0 - window
    window_end = pos0 + window

    base_at_target = None
    window_items = []

    ref_pos = ref_start
    query_pos = 0
    prev_ref_pos = None
    read_len = len(read_bases)

    # pysam CIGAR operation codes:
    # 0 M, 1 I, 2 D, 3 N, 4 S, 5 H, 6 P, 7 =, 8 X
    for op, length in cigartuples:
        if length <= 0:
            continue

        if op == 0 or op == 7 or op == 8:
            block_ref_start = ref_pos
            block_ref_end = ref_pos + length - 1

            if block_ref_end >= window_start and block_ref_start <= window_end:
                overlap_start = max(block_ref_start, window_start)
                overlap_end = min(block_ref_end, window_end)

                q = query_pos + (overlap_start - block_ref_start)

                for r in range(overlap_start, overlap_end + 1):
                    window_idx = r - window_start
                    window_items.append((window_idx, q, r))

                    if r == pos0 and q < read_len:
                        base_at_target = read_bases[q].upper()

                    q += 1

            prev_ref_pos = block_ref_end
            ref_pos += length
            query_pos += length

            if ref_pos > window_end + 1:
                break

        elif op == 1:
            if prev_ref_pos is not None:
                if window_start <= prev_ref_pos <= window_end:
                    window_idx = prev_ref_pos - window_start

                    for _ in range(length):
                        window_items.append((window_idx, "INS", prev_ref_pos))

                if prev_ref_pos == pos0:
                    base_at_target = "INS"

            query_pos += length

        elif op == 2 or op == 3:
            block_ref_start = ref_pos
            block_ref_end = ref_pos + length - 1

            if block_ref_end >= window_start and block_ref_start <= window_end:
                overlap_start = max(block_ref_start, window_start)
                overlap_end = min(block_ref_end, window_end)

                for r in range(overlap_start, overlap_end + 1):
                    window_idx = r - window_start
                    window_items.append((window_idx, None, r))

                    if r == pos0:
                        base_at_target = "DEL"

            prev_ref_pos = block_ref_end
            ref_pos += length

            if ref_pos > window_end + 1:
                break

        elif op == 4:
            query_pos += length

        elif op == 5 or op == 6:
            continue

        else:
            continue

    return base_at_target, window_items


def extract_features_for_snp_from_reads(
    reads,
    chrom,
    pos,
    ref,
    alt,
):
    """
    Extract features for one SNP from a pre-selected list of reads.

    This is the core feature function used by the chunk-based extractor.
    It does not fetch reads from BAM.
    """
    ref = ref.upper()
    alt = alt.upper()

    reference_sequence = get_reference_sequence(chrom, pos)
    one_hot_encoded_sequence = one_hot_encode_sequence(reference_sequence)

    acc = {
        "total": init_feature_accumulators(),
        "ref_forward": init_feature_accumulators(),
        "ref_reverse": init_feature_accumulators(),
        "alt_forward": init_feature_accumulators(),
        "alt_reverse": init_feature_accumulators(),
    }

    for read in reads:
        # Safety check. Most filtering is already done at chunk-fetch level.
        if not read_passes_basic_filter(read, min_mapq=20):
            continue

        base_at_target, window_items = extract_read_target_and_window(read, pos)

        if window_items is None:
            continue

        group = None

        if base_at_target == ref:
            group = "ref_reverse" if read.is_reverse else "ref_forward"
        elif base_at_target == alt:
            group = "alt_reverse" if read.is_reverse else "alt_forward"

        if group is None:
            continue

        add_read_to_accumulator(acc[group], read, window_items, pos)
        add_read_to_accumulator(acc["total"], read, window_items, pos)

    total_base_fraction, total_coverage = finalize_base_fraction_and_coverage(acc["total"])

    target_position_coverage = np.array(
        [min(float(total_coverage[30]) / 100.0, 1.0)],
        dtype=np.float32
    )

    ref_base_fractions_forward, ref_coverage_forward = finalize_base_fraction_and_coverage(
        acc["ref_forward"]
    )
    ref_base_fractions_reverse, ref_coverage_reverse = finalize_base_fraction_and_coverage(
        acc["ref_reverse"]
    )
    alt_base_fractions_forward, alt_coverage_forward = finalize_base_fraction_and_coverage(
        acc["alt_forward"]
    )
    alt_base_fractions_reverse, alt_coverage_reverse = finalize_base_fraction_and_coverage(
        acc["alt_reverse"]
    )

    scale_coverage = max(
        float(np.max(ref_coverage_forward)),
        float(np.max(ref_coverage_reverse)),
        float(np.max(alt_coverage_forward)),
        float(np.max(alt_coverage_reverse)),
    )

    if scale_coverage > 0:
        ref_coverage_forward = ref_coverage_forward / scale_coverage
        ref_coverage_reverse = ref_coverage_reverse / scale_coverage
        alt_coverage_forward = alt_coverage_forward / scale_coverage
        alt_coverage_reverse = alt_coverage_reverse / scale_coverage

    ref_mapq_forward = finalize_average_mapq(acc["ref_forward"])
    ref_mapq_reverse = finalize_average_mapq(acc["ref_reverse"])
    alt_mapq_forward = finalize_average_mapq(acc["alt_forward"])
    alt_mapq_reverse = finalize_average_mapq(acc["alt_reverse"])

    ref_baseq_forward = finalize_average_baseq(acc["ref_forward"])
    ref_baseq_reverse = finalize_average_baseq(acc["ref_reverse"])
    alt_baseq_forward = finalize_average_baseq(acc["alt_forward"])
    alt_baseq_reverse = finalize_average_baseq(acc["alt_reverse"])

    ref_mismatch_rate_forward = finalize_average_mismatch_rate(acc["ref_forward"])
    ref_mismatch_rate_reverse = finalize_average_mismatch_rate(acc["ref_reverse"])
    alt_mismatch_rate_forward = finalize_average_mismatch_rate(acc["alt_forward"])
    alt_mismatch_rate_reverse = finalize_average_mismatch_rate(acc["alt_reverse"])

    ref_read_length_forward = finalize_average_read_length(acc["ref_forward"])
    ref_read_length_reverse = finalize_average_read_length(acc["ref_reverse"])
    alt_read_length_forward = finalize_average_read_length(acc["alt_forward"])
    alt_read_length_reverse = finalize_average_read_length(acc["alt_reverse"])

    ref_target_relative_position_forward = finalize_average_target_relative_position(acc["ref_forward"])
    ref_target_relative_position_reverse = finalize_average_target_relative_position(acc["ref_reverse"])
    alt_target_relative_position_forward = finalize_average_target_relative_position(acc["alt_forward"])
    alt_target_relative_position_reverse = finalize_average_target_relative_position(acc["alt_reverse"])

    ref_target_coverage_forward = finalize_target_coverage_count(acc["ref_forward"])
    ref_target_coverage_reverse = finalize_target_coverage_count(acc["ref_reverse"])
    alt_target_coverage_forward = finalize_target_coverage_count(acc["alt_forward"])
    alt_target_coverage_reverse = finalize_target_coverage_count(acc["alt_reverse"])

    return (
        one_hot_encoded_sequence,
        total_base_fraction,
        target_position_coverage,

        ref_base_fractions_forward,
        ref_base_fractions_reverse,
        alt_base_fractions_forward,
        alt_base_fractions_reverse,

        ref_coverage_forward,
        ref_coverage_reverse,
        alt_coverage_forward,
        alt_coverage_reverse,

        ref_mapq_forward,
        ref_mapq_reverse,
        alt_mapq_forward,
        alt_mapq_reverse,

        ref_baseq_forward,
        ref_baseq_reverse,
        alt_baseq_forward,
        alt_baseq_reverse,

        ref_mismatch_rate_forward,
        ref_mismatch_rate_reverse,
        alt_mismatch_rate_forward,
        alt_mismatch_rate_reverse,

        ref_read_length_forward,
        ref_read_length_reverse,
        alt_read_length_forward,
        alt_read_length_reverse,

        ref_target_relative_position_forward,
        ref_target_relative_position_reverse,
        alt_target_relative_position_forward,
        alt_target_relative_position_reverse,

        ref_target_coverage_forward,
        ref_target_coverage_reverse,
        alt_target_coverage_forward,
        alt_target_coverage_reverse,
    )


def extract_features_for_chunk_wrapper(args):
    """
    Process one genomic chunk.

    For all variants in a chunk:
      1. Fetch reads once across the chunk region.
      2. Filter reads once.
      3. Assign each read to every candidate position it overlaps.
      4. Apply per-position read_depth cap.
      5. Extract features for each variant from the assigned reads.
    """
    (
        chunk_id,
        bam_file,
        reference_file,
        chrom,
        fetch_start,
        fetch_end,
        variants,
        read_depth,
    ) = args

    global GLOBAL_BAM

    variants = sorted(variants, key=lambda x: x["pos"])
    pos0_list = [int(v["pos"]) - 1 for v in variants]
    reads_by_variant = [[] for _ in variants]

    n_variants = len(variants)
    full_variant_count = 0
    is_full = [False] * n_variants

    for read in GLOBAL_BAM.fetch(chrom, fetch_start, fetch_end):
        if not read_passes_basic_filter(read, min_mapq=20):
            continue

        read_start = read.reference_start
        read_end = read.reference_end

        if read_start is None or read_end is None:
            continue

        # Find candidate positions with pos0 >= read_start.
        j = bisect_left(pos0_list, read_start)

        # Add read to all candidate positions with pos0 < read_end.
        while j < n_variants and pos0_list[j] < read_end:
            if read_depth <= 0 or len(reads_by_variant[j]) < read_depth:
                reads_by_variant[j].append(read)

                if read_depth > 0 and len(reads_by_variant[j]) >= read_depth and not is_full[j]:
                    is_full[j] = True
                    full_variant_count += 1

            j += 1

        # If every candidate position in this chunk already has enough reads,
        # stop reading more BAM records for this chunk.
        if read_depth > 0 and full_variant_count >= n_variants:
            break

    chunk_results = []

    for v, reads in zip(variants, reads_by_variant):
        features = extract_features_for_snp_from_reads(
            reads=reads,
            chrom=v["chrom"],
            pos=int(v["pos"]),
            ref=v["ref"],
            alt=v["alt"],
        )

        chunk_results.append((
            v["index"],
            v["chrom"],
            int(v["pos"]),
            v["ref"],
            v["alt"],
            features,
        ))

    return chunk_id, chunk_results


def feature_tuple_to_dict(features):
    return {
        "one_hot_encoded_sequence": features[0],
        "total_base_fraction": features[1],
        "target_position_coverage": features[2],

        "ref_base_fractions_forward": features[3],
        "ref_base_fractions_reverse": features[4],
        "alt_base_fractions_forward": features[5],
        "alt_base_fractions_reverse": features[6],

        "ref_coverage_forward": features[7],
        "ref_coverage_reverse": features[8],
        "alt_coverage_forward": features[9],
        "alt_coverage_reverse": features[10],

        "ref_mapq_forward": features[11],
        "ref_mapq_reverse": features[12],
        "alt_mapq_forward": features[13],
        "alt_mapq_reverse": features[14],

        "ref_baseq_forward": features[15],
        "ref_baseq_reverse": features[16],
        "alt_baseq_forward": features[17],
        "alt_baseq_reverse": features[18],

        "ref_mismatch_rate_forward": features[19],
        "ref_mismatch_rate_reverse": features[20],
        "alt_mismatch_rate_forward": features[21],
        "alt_mismatch_rate_reverse": features[22],

        "ref_read_length_forward": features[23],
        "ref_read_length_reverse": features[24],
        "alt_read_length_forward": features[25],
        "alt_read_length_reverse": features[26],

        "ref_target_relative_position_forward": features[27],
        "ref_target_relative_position_reverse": features[28],
        "alt_target_relative_position_forward": features[29],
        "alt_target_relative_position_reverse": features[30],

        "ref_target_coverage_forward": features[31],
        "ref_target_coverage_reverse": features[32],
        "alt_target_coverage_forward": features[33],
        "alt_target_coverage_reverse": features[34],
    }


def prepare_feature_arrays(data):
    arrays = {}
    for name in FEATURE_ORDER:
        arrays[name] = np.array([d[name] for d in data], dtype=np.float32)
    return arrays


class H5ShardWriter:
    def __init__(self, output_dir, shard_size=100000, compression="lzf"):
        self.output_dir = output_dir
        self.shard_size = shard_size
        self.compression = compression
        self.buffer = []
        self.shard_id = 0
        self.total_examples = 0
        self.manifest_rows = []

        os.makedirs(output_dir, exist_ok=True)

    def add(self, record):
        self.buffer.append(record)

        if len(self.buffer) >= self.shard_size:
            self.flush()

    def flush(self):
        if len(self.buffer) == 0:
            return

        shard_name = f"shard_{self.shard_id:06d}.h5"
        shard_path = os.path.join(self.output_dir, shard_name)

        feature_data = [r["features"] for r in self.buffer]
        arrays = prepare_feature_arrays(feature_data)

        chrom = np.array([r["chrom"] for r in self.buffer], dtype=object)
        pos = np.array([r["pos"] for r in self.buffer], dtype=np.int64)
        ref = np.array([r["ref"] for r in self.buffer], dtype=object)
        alt = np.array([r["alt"] for r in self.buffer], dtype=object)
        key = np.array(
            [f'{r["chrom"]}:{r["pos"]}:{r["ref"]}:{r["alt"]}' for r in self.buffer],
            dtype=object
        )

        str_dtype = h5py.string_dtype(encoding="utf-8")

        with h5py.File(shard_path, "w") as h5:
            h5.attrs["num_examples"] = len(self.buffer)
            h5.attrs["feature_order"] = ",".join(FEATURE_ORDER)
            h5.attrs["base_fraction_row_order"] = "A,C,G,T,DEL,INS"

            for name in FEATURE_ORDER:
                h5.create_dataset(
                    name,
                    data=arrays[name],
                    compression=self.compression
                )

            h5.create_dataset("chrom", data=chrom.astype(str_dtype), dtype=str_dtype)
            h5.create_dataset("pos", data=pos)
            h5.create_dataset("ref", data=ref.astype(str_dtype), dtype=str_dtype)
            h5.create_dataset("alt", data=alt.astype(str_dtype), dtype=str_dtype)
            h5.create_dataset("key", data=key.astype(str_dtype), dtype=str_dtype)

        for local_idx, r in enumerate(self.buffer):
            self.manifest_rows.append([
                shard_name,
                local_idx,
                r["chrom"],
                r["pos"],
                r["ref"],
                r["alt"],
                f'{r["chrom"]}:{r["pos"]}:{r["ref"]}:{r["alt"]}',
            ])

        print(f"Wrote {shard_path} with {len(self.buffer)} examples")

        self.total_examples += len(self.buffer)
        self.buffer = []
        self.shard_id += 1

    def close(self):
        self.flush()

        manifest_path = os.path.join(self.output_dir, "manifest.tsv")

        with open(manifest_path, "w") as f:
            f.write("shard\tlocal_idx\tchrom\tpos\tref\talt\tkey\n")
            for row in self.manifest_rows:
                f.write("\t".join(map(str, row)) + "\n")

        print(f"Wrote manifest: {manifest_path}")
        print(f"Total examples: {self.total_examples}")


def extract_VCF_feature(
    vcf_file,
    bam_file,
    reference_file,
    num_processes,
    read_depth,
    shard_writer,
    chunk_bp=1000,
):
    chunk_tasks, total_variants = make_chunk_tasks(
        vcf_file=vcf_file,
        bam_file=bam_file,
        reference_file=reference_file,
        read_depth=read_depth,
        chunk_bp=chunk_bp,
    )

    if len(chunk_tasks) == 0:
        print(f"No variants to process in {vcf_file}")
        return

    print(
        f"Chunk-based extraction: {total_variants} variants in "
        f"{len(chunk_tasks)} chunks; chunk_bp={chunk_bp}; read_depth={read_depth}"
    )

    widgets = [
        "Extract SNV features: ",
        progressbar.Percentage(), " ",
        progressbar.Bar(marker="=", left="[", right="]"), " ",
        progressbar.ETA()
    ]
    bar = progressbar.ProgressBar(widgets=widgets, maxval=total_variants).start()

    count = 0

    with Pool(
        processes=num_processes,
        initializer=init_worker,
        initargs=(bam_file, reference_file)
    ) as pool:

        for _chunk_id, chunk_results in pool.imap(
            extract_features_for_chunk_wrapper,
            chunk_tasks,
            chunksize=1,
        ):
            # chunk_results are already sorted by position within chunk.
            # chunk_tasks are processed in sorted chunk order by imap.
            for (
                _index,
                chrom,
                pos,
                ref,
                alt,
                features,
            ) in chunk_results:
                feature_dict = feature_tuple_to_dict(features)

                shard_writer.add({
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "features": feature_dict,
                })

                count += 1
                bar.update(count)

    bar.finish()


def run_snp_feature_extraction_multiple(args):
    filtered_variant = filter_vcf(
        args.variant,
        args.filtered_variant,
        args.min_vaf,
        args.min_total_coverage,
        args.ALT,
    )

    shard_writer = H5ShardWriter(
        output_dir=args.output,
        shard_size=args.shard_size,
        compression=args.compression,
    )

    if args.all_v:
        filtered_variant.to_csv(args.filtered_variant, sep="\t", index=False)

        extract_VCF_feature(
            args.filtered_variant,
            args.bam,
            args.ref,
            args.threads,
            args.depth,
            shard_writer,
            chunk_bp=args.chunk_bp,
        )

    else:
        filtered_variant.to_csv(args.filtered_variant, sep="\t", index=False)

        for i in range(1, 23):
            each_f_var_file = args.filtered_variant + f".chr{i}"

            each_filtered_variant = filtered_variant.loc[
                filtered_variant.loc[:, "chrom"] == f"chr{i}",
                :
            ]

            if each_filtered_variant.shape[0] == 0:
                continue

            print(f"Start processing chr{i}")

            each_filtered_variant.to_csv(each_f_var_file, sep="\t", index=False)

            extract_VCF_feature(
                each_f_var_file,
                args.bam,
                args.ref,
                args.threads,
                args.depth,
                shard_writer,
                chunk_bp=args.chunk_bp,
            )

    shard_writer.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract simplified SNP features from VCF/BAM and write HDF5 feature shards."
    )

    parser.add_argument("-v", "--variant", required=True, help="Path to input SNV file")
    parser.add_argument("--bam", required=True, help="Path to input BAM file")
    parser.add_argument("--ref", required=True, help="Path to reference genome FASTA file")
    parser.add_argument("--output", required=True, help="Output directory for HDF5 shards")
    parser.add_argument("--all_v", help="Kept for compatibility; not used in simplified feature extraction")
    parser.add_argument("--filtered_variant", required=True, help="Path to output filtered SNV file")

    parser.add_argument("-t", "--threads", type=int, default=24)
    parser.add_argument("-d", "--depth", type=int, default=1000)

    parser.add_argument("--min-vaf", type=float, default=0.05)
    parser.add_argument("--min-total-coverage", type=int, default=3)
    parser.add_argument("--ALT", type=int, default=3)

    parser.add_argument(
        "--chunk-bp",
        type=int,
        default=1000,
        help="Genomic chunk size for chunk-based BAM fetching. Default: 1000."
    )

    parser.add_argument(
        "--shard_size",
        type=int,
        default=100000,
        help="Number of examples per HDF5 shard. Default: 100000"
    )

    parser.add_argument(
        "--compression",
        default="lzf",
        choices=["lzf", "gzip", None],
        help="HDF5 compression. lzf is faster; gzip is smaller."
    )

    args = parser.parse_args()

    run_snp_feature_extraction_multiple(args)
