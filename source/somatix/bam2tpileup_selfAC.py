import sys
import argparse
import subprocess
import re
import pysam
from multiprocessing import Pool
import progressbar


def parse_allele_counter_output(input_stream):
    """
    Expected updated allele_counter columns:
    chrom pos ref alt ref_count alt_count total_coverage vaf A C G T INS DEL

    Output columns:
    chrom pos ref alt ref_reads alt_reads total_coverage vaf A C G T INS DEL alt_tie

    ALT selection:
    1. Choose the non-reference allele with maximum read count.
    2. If tied, choose by order: INS, DEL, A, C, G, T.

    alt_tie:
    All non-reference alleles with the same read count as the selected ALT,
    ordered as  A, C, G, T, INS, DEL.
    """
    results = []
    alt_order = ["A", "T", "C", "G",  "INS", "DEL"]

    for line in input_stream:
        line = line.strip()
        if not line or line.startswith("chrom\t"):
            continue

        cols = line.split("\t")
        if len(cols) < 14:
            continue

        chrom = cols[0]
        pos = int(cols[1])
        ref_base = cols[2].upper()

        ref_count = int(cols[4])
        total_coverage = int(cols[6])

        A = int(cols[8])
        C = int(cols[9])
        G = int(cols[10])
        T = int(cols[11])
        INS = int(cols[12])
        DEL = int(cols[13])

        allele_counts = {
            "INS": INS,
            "DEL": DEL,
            "A": A,
            "C": C,
            "G": G,
            "T": T,
        }

        candidate_alleles = [
            allele for allele in alt_order
            if allele != ref_base
        ]

        max_alt_count = max(allele_counts[allele] for allele in candidate_alleles)

        if max_alt_count <= 0:
            continue

        alt = next(
            allele for allele in candidate_alleles
            if allele_counts[allele] == max_alt_count
        )

        alt_count = max_alt_count
        vaf = alt_count / total_coverage if total_coverage > 0 else 0

        alt_tie = [
            allele for allele in candidate_alleles
            if allele_counts[allele] == alt_count and alt_count > 0
        ]

        results.append("\t".join(str(x) for x in [
            chrom, pos, ref_base, alt,
            ref_count, alt_count,
            total_coverage, vaf,
            A, C, G, T, INS, DEL,
            ",".join(alt_tie)
        ]))

    return results


def run_allele_counter_wrapper(args):
    return run_allele_counter(*args)


def run_allele_counter(
    region,
    bam_path,
    ref_path,
    allele_counter_path,
    min_alt,
    min_vaf,
    min_total_coverage,
    min_mapq,
    min_baseq,
    max_depth,
    excl_flags,
):
    cmd = [
        allele_counter_path,
        "--bam", bam_path,
        "--ref", ref_path,
        "--region", region,
        "--min-mapq", str(min_mapq),
        "--min-baseq", str(min_baseq),
        "--min-alt", str(min_alt),
        "--min-total-coverage", str(min_total_coverage),
        "--min-vaf", str(min_vaf),
        "--max-depth", str(max_depth),
        "--excl-flags", str(excl_flags),
    ]

    process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    stdout, stderr = process.communicate()

    if process.returncode != 0:
        raise RuntimeError(
            f"allele_counter failed for {region}\n"
            f"Command: {' '.join(cmd)}\n"
            f"stderr:\n{stderr}"
        )

    return parse_allele_counter_output(stdout.splitlines())


def split_regions(chrom, start_pos, end_pos, block_size=1_000_000):
    regions = []
    for start in range(start_pos, end_pos + 1, block_size):
        end = min(start + block_size - 1, end_pos)
        regions.append(f"{chrom}:{start}-{end}")
    return regions


def process_bam(
    bam_path,
    ref_path,
    allele_counter_path,
    region=None,
    min_alt=2,
    min_vaf=0.05,
    min_total_coverage=10,
    threads=4,
    block_size=1_000_000,
    min_mapq=20,
    min_baseq=10,
    max_depth=5000,
    excl_flags=2316,
):
    bam = pysam.AlignmentFile(bam_path, "rb")
    chrom_lengths = {ref["SN"]: ref["LN"] for ref in bam.header["SQ"]}
    bam.close()

    headline = "\t".join([
        "chrom", "pos", "ref", "alt",
        "ref_reads", "alt_reads",
        "total_coverage", "vaf",
        "A", "C", "G", "T", "INS", "DEL",
        "alt_tie",
    ])
    sys.stdout.write(headline + "\n")

    regions = []

    if region:
        match = re.match(r"([^:]+)(?::(\d+)-(\d+))?$", region)
        if not match:
            raise ValueError("Invalid region format. Use 'chrom:start-end' or 'chrom'.")

        chrom = match.group(1)
        if chrom not in chrom_lengths:
            raise ValueError(f"Chromosome {chrom} not found in BAM header.")

        start = int(match.group(2)) if match.group(2) else 1
        end = int(match.group(3)) if match.group(3) else chrom_lengths[chrom]

        if start < 1 or start > end:
            raise ValueError("Invalid region start/end.")
        if end > chrom_lengths[chrom]:
            end = chrom_lengths[chrom]

        regions = split_regions(chrom, start, end, block_size=block_size)

    else:
        target_chrom = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
        for chrom in target_chrom:
            if chrom in chrom_lengths:
                regions.extend(
                    split_regions(chrom, 1, chrom_lengths[chrom], block_size=block_size)
                )

    widgets = [
        "Detect variants: ",
        progressbar.Percentage(), " ",
        progressbar.Bar(marker="=", left="[", right="]"), " ",
        progressbar.ETA(),
    ]
    bar = progressbar.ProgressBar(widgets=widgets, maxval=len(regions)).start()

    tasks = [
        (
            r,
            bam_path,
            ref_path,
            allele_counter_path,
            min_alt,
            min_vaf,
            min_total_coverage,
            min_mapq,
            min_baseq,
            max_depth,
            excl_flags,
        )
        for r in regions
    ]

    count = 0
    with Pool(processes=threads) as pool:
        for result in pool.imap(run_allele_counter_wrapper, tasks):
            for line in result:
                sys.stdout.write(line + "\n")
            count += 1
            bar.update(count)

    bar.finish()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract candidate variant sites using allele_counter backend."
    )

    parser.add_argument("--bam", required=True)
    parser.add_argument("--ref", required=True)
    parser.add_argument("--region", required=False)

    parser.add_argument(
        "--allele-counter",
        default="allele_counter",
        help="Path to compiled allele_counter binary",
    )

    parser.add_argument(
        "--min-alt",
        type=int,
        default=2,
        help="Minimum reported ALT read count",
    )

    parser.add_argument(
        "--min-vaf",
        type=float,
        default=0.05,
        help="Minimum ALT VAF = alt_reads / total_coverage. Default: 0.05",
    )

    parser.add_argument(
        "--min-total-coverage",
        type=int,
        default=10,
        help="Minimum total coverage using A+C+G+T+INS+DEL. Default: 10",
    )

    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--block-size", type=int, default=1_000_000)

    parser.add_argument("--min-mapq", "--min-MQ", type=int, default=20)
    parser.add_argument("--min-baseq", "--min-BQ", type=int, default=10)
    parser.add_argument("--max-depth", type=int, default=5000)
    parser.add_argument("--excl-flags", type=int, default=2316)

    args = parser.parse_args()

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
