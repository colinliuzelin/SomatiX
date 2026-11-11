#!/usr/bin/env python3
"""Exercise the container with synthetic reads only; no benchmark downloads."""

import argparse
import importlib
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


MODEL_NAMES = (
    "somatix_ont_multicancer.pth",
    "somatix_ont_hcc1395.pth",
    "somatix_pacbio_multicancer.pth",
    "somatix_pacbio_hcc1395.pth",
)


def run(command):
    print("[TEST]", " ".join(map(str, command)), flush=True)
    subprocess.run(list(map(str, command)), check=True, timeout=180)


def write_bam(path, reference, alt_count):
    import pysam

    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": "chr1", "LN": len(reference)}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for index in range(40):
            start = 100 + index
            sequence = list(reference[start:start + 700])
            if index < alt_count:
                sequence[500 - start] = "T"
            read = pysam.AlignedSegment(bam.header)
            read.query_name = f"read_{index}"
            read.query_sequence = "".join(sequence)
            read.flag = 16 if index % 2 else 0
            read.reference_id = 0
            read.reference_start = start
            read.mapping_quality = 60
            read.cigar = [(0, len(sequence))]
            read.query_qualities = pysam.qualitystring_to_array("I" * len(sequence))
            read.set_tag("NM", int(index < alt_count))
            bam.write(read)
    pysam.index(str(path))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path("/opt/SomatiX"))
    parser.add_argument("--somatix", nargs="+", default=["somatix"])
    parser.add_argument("--allele-counter", default="allele_counter")
    parser.add_argument("--require-cpu", action="store_true")
    args = parser.parse_args()
    sys.path.insert(0, str(args.repo / "source"))
    for name in ("numpy", "pandas", "h5py", "pysam", "progressbar", "torch"):
        importlib.import_module(name)
    import numpy as np
    import pandas as pd
    import pysam
    import torch
    from somatix.predict_DNN_somatic import DNASeqCoverageModel_somatic, predict_one_chrom

    torch.set_num_threads(1)
    if args.require_cpu and torch.version.cuda is not None:
        raise RuntimeError("Expected CPU-only PyTorch, but a CUDA build was installed")
    counter = shutil.which(args.allele_counter)
    if counter is None:
        raise FileNotFoundError(args.allele_counter)
    libraries = subprocess.run(["ldd", counter], capture_output=True, text=True, check=True)
    if "not found" in libraries.stdout + libraries.stderr:
        raise RuntimeError(libraries.stdout + libraries.stderr)
    run([counter, "--help"])
    run(args.somatix + ["--version"])
    run(args.somatix + ["call", "--help"])

    with tempfile.TemporaryDirectory(prefix="somatix-smoke-") as directory:
        root = Path(directory)
        reference = "ACGT" * 300
        fasta = root / "reference.fa"
        fasta.write_text(">chr1\n" + reference + "\n", encoding="ascii")
        pysam.faidx(str(fasta))
        tumor, normal = root / "tumor.bam", root / "normal.bam"
        write_bam(tumor, reference, alt_count=12)
        write_bam(normal, reference, alt_count=0)
        output = root / "output"
        run(args.somatix + [
            "call", "--bam-case", tumor, "--bam-control", normal,
            "--ref", fasta, "--model", args.repo / "model" / MODEL_NAMES[0],
            "--outdir", output, "--region", "chr1", "--threads", "1",
            "--device", "cpu", "--allele-counter", counter,
        ])
        predictions = pd.read_csv(output / "somatix_predict.txt", sep="\t")
        assert len(predictions) == 1, "Expected exactly one synthetic candidate"
        assert int(predictions.iloc[0]["pos"]) == 501
        assert predictions.iloc[0]["ref"] == "A" and predictions.iloc[0]["alt"] == "T"
        for filename in ("somatix_predict.vcf", "somatix_predict.somatic.vcf",
                         "somatix_predict.germline.vcf"):
            with pysam.VariantFile(str(output / filename)) as vcf:
                list(vcf)

        for name in MODEL_NAMES:
            state = torch.load(args.repo / "model" / name, map_location="cpu", weights_only=True)
            if isinstance(state, dict) and "state_dict" in state:
                state = state["state_dict"]
            model = DNASeqCoverageModel_somatic()
            model.load_state_dict(state, strict=True)
            model.eval()
            result = predict_one_chrom(
                str(output / "case_features.chr1"),
                str(output / "control_features.chr1"),
                str(output / "case_candidates.filtered.chr1"),
                model, torch.device("cpu"),
            )
            probabilities = result[[f"class_{i}" for i in range(4)]].to_numpy()
            assert probabilities.shape == (1, 4), name
            assert np.isfinite(probabilities).all(), name
            assert np.allclose(probabilities.sum(axis=1), 1.0), name
            print(f"[PASS] Checkpoint and forward pass: {name}", flush=True)
    print("[PASS] Libraries, CLI, all checkpoints and full calling workflow", flush=True)


if __name__ == "__main__":
    for variable in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
        os.environ.setdefault(variable, "1")
    main()
