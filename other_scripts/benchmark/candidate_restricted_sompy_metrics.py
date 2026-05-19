#!/usr/bin/env python3
"""Calculate candidate-restricted som.py SNV metrics.

This script keeps the benchmarking logic tied to som.py: TP, FP and FN are
taken from the som.py feature-table tags after restricting rows to SomatiX
candidate sites, using chrom:pos as the site key.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


SNV_BASES = {"A", "C", "G", "T"}


def normalize_chrom(chrom: str) -> str:
    chrom = str(chrom).strip()
    if not chrom:
        return chrom
    if chrom.startswith("chr"):
        return chrom
    return f"chr{chrom}"


def site_key(chrom: str, pos: str | int) -> str:
    return f"{normalize_chrom(chrom)}:{int(pos)}"


def is_snv(ref: str, alt: str) -> bool:
    return str(ref).upper() in SNV_BASES and str(alt).upper() in SNV_BASES


def as_float(value: str, default: float = math.nan) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def read_candidate_sites(
    path: Path,
    min_alt: int,
    min_vaf: float,
    chrom: str | None,
) -> set[str]:
    keys: set[str] = set()
    wanted_chrom = normalize_chrom(chrom) if chrom else None

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"chrom", "pos", "ref", "alt", "alt_reads", "vaf"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"Candidate file is missing required columns: {', '.join(sorted(missing))}"
            )

        for row in reader:
            row_chrom = normalize_chrom(row["chrom"])
            if wanted_chrom and row_chrom != wanted_chrom:
                continue
            if not is_snv(row["ref"], row["alt"]):
                continue
            if as_float(row["alt_reads"]) < min_alt:
                continue
            if as_float(row["vaf"]) < min_vaf:
                continue
            keys.add(site_key(row_chrom, row["pos"]))

    return keys


def sompy_row_is_snv(row: dict[str, str]) -> bool:
    tag = row.get("tag", "")
    if tag == "FN":
        return is_snv(row.get("REF.truth", ""), row.get("ALT.truth", ""))
    return is_snv(row.get("REF", ""), row.get("ALT", ""))


def summarize_sompy_features(
    path: Path,
    candidate_keys: set[str],
    chrom: str | None,
    filtered_output: Path | None,
) -> dict[str, int | float | str]:
    wanted_chrom = normalize_chrom(chrom) if chrom else None
    counts = {"TP": 0, "FP": 0, "FN": 0, "UNK": 0, "AMBI": 0, "AMBIG": 0}
    sompy_rows_total = 0
    sompy_snv_rows_total = 0
    sompy_rows_in_candidates = 0

    out_handle = None
    writer = None
    try:
        if filtered_output:
            filtered_output.parent.mkdir(parents=True, exist_ok=True)
            out_handle = filtered_output.open("w", newline="")

        with path.open(newline="") as handle:
            reader = csv.DictReader(handle)
            required = {"CHROM", "POS", "tag"}
            missing = required.difference(reader.fieldnames or [])
            if missing:
                raise ValueError(
                    f"som.py feature table is missing required columns: {', '.join(sorted(missing))}"
                )

            if out_handle:
                writer = csv.DictWriter(out_handle, fieldnames=reader.fieldnames)
                writer.writeheader()

            for row in reader:
                sompy_rows_total += 1
                row_chrom = normalize_chrom(row["CHROM"])
                if wanted_chrom and row_chrom != wanted_chrom:
                    continue
                if not sompy_row_is_snv(row):
                    continue
                sompy_snv_rows_total += 1
                if site_key(row_chrom, row["POS"]) not in candidate_keys:
                    continue

                sompy_rows_in_candidates += 1
                tag = row.get("tag", "")
                if tag in counts:
                    counts[tag] += 1
                if writer:
                    writer.writerow(row)
    finally:
        if out_handle:
            out_handle.close()

    tp = counts["TP"]
    fp = counts["FP"]
    fn = counts["FN"]
    precision = tp / (tp + fp) if (tp + fp) > 0 else math.nan
    recall = tp / (tp + fn) if (tp + fn) > 0 else math.nan
    f1 = (
        2 * precision * recall / (precision + recall)
        if math.isfinite(precision) and math.isfinite(recall) and (precision + recall) > 0
        else math.nan
    )

    return {
        "candidate_sites": len(candidate_keys),
        "sompy_rows_total": sompy_rows_total,
        "sompy_snv_rows_total": sompy_snv_rows_total,
        "sompy_rows_in_candidates": sompy_rows_in_candidates,
        "total_truth": tp + fn,
        "total_query": tp + fp,
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "unk": counts["UNK"],
        "ambi": counts["AMBI"] + counts["AMBIG"],
        "precision": precision,
        "recall": recall,
        "f1": f1,
    }


def write_metrics(path: Path, metrics: dict[str, int | float | str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(metrics))
        writer.writeheader()
        writer.writerow(metrics)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Calculate candidate-restricted som.py SNV metrics."
    )
    parser.add_argument("--sompy-features", required=True, type=Path)
    parser.add_argument("--candidates", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--filtered-features-output", type=Path)
    parser.add_argument("--sample-id", default="")
    parser.add_argument("--platform", default="")
    parser.add_argument("--tool", default="")
    parser.add_argument("--data-type", default="")
    parser.add_argument("--chrom", default=None)
    parser.add_argument("--min-alt", type=int, default=3)
    parser.add_argument("--min-vaf", type=float, default=0.05)
    args = parser.parse_args()

    candidate_keys = read_candidate_sites(
        args.candidates,
        min_alt=args.min_alt,
        min_vaf=args.min_vaf,
        chrom=args.chrom,
    )
    summary = summarize_sompy_features(
        args.sompy_features,
        candidate_keys=candidate_keys,
        chrom=args.chrom,
        filtered_output=args.filtered_features_output,
    )

    metrics = {
        "sample_id": args.sample_id,
        "platform": args.platform,
        "tool": args.tool,
        "data_type": args.data_type,
        "chrom": normalize_chrom(args.chrom) if args.chrom else "",
        "min_alt": args.min_alt,
        "min_vaf": args.min_vaf,
        "candidate_file": str(args.candidates),
        "sompy_features": str(args.sompy_features),
        **summary,
    }
    write_metrics(args.output, metrics)

    print(
        "candidate_sites={candidate_sites} tp={tp} fp={fp} fn={fn} "
        "precision={precision:.6g} recall={recall:.6g} f1={f1:.6g}".format(
            **metrics
        )
    )
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
