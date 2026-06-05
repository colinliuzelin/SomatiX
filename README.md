# SomatiX

SomatiX is a deep-learning variant caller for accurate and fast somatic SNV detection from **long-read DNA sequencing (tumor-normal pair samples)**. It supports both **Oxford Nanopore (ONT)** and **Pacific Biosciences (PacBio)** platforms, with pretrained models for each technology type.

## Contents

- [Installation](#installation)
- [Case Study](#case-study)
- [Quick Start and Overview](#quick-start-and-overview)
- [Model Architecture](#model-architecture)
- [Pretrained Models](#pretrained-models)
- [Basic Usage](#basic-usage)
- [Detailed CLI Options](#detailed-cli-options)
  - [`somatix candidates`](#somatix-candidates)
  - [`somatix features`](#somatix-features)
  - [`somatix predict`](#somatix-predict)
  - [`somatix perm`](#somatix-perm)
  - [`somatix call`](#somatix-call)
  - [`somatix clean`](#somatix-clean)
- [Outputs](#outputs)
- [Resource Notes](#resource-notes)
- [License](#license)
- [Contact](#contact)

---

## Installation

Supported OS: Linux.

### Option 1: Conda
First clone the repository, or download it from GitHub and enter the repository
root:

```bash
git clone git@github.com:colinliuzelin/SomatiX.git
cd SomatiX
```

Ensure you have **Conda** installed, then create and activate an environment
from the included `environment.yml`:

```bash
conda env create -f environment.yml -n somatix
conda activate somatix

# IMPORTANT: `environment.yml` does NOT install PyTorch by default. Install
# PyTorch separately to match your system CUDA runtime (or use a CPU-only
# build). Example commands:

```bash
# CPU-only (conda, official PyTorch channel)
conda install -y pytorch torchvision torchaudio cpuonly -c pytorch

# CUDA 12.1 (example; change to the CUDA version you have)
conda install -y pytorch torchvision torchaudio pytorch-cuda=12.1 -c pytorch -c nvidia

# Or use pip with the appropriate wheel index (example for CUDA 11.8):
pip install --index-url https://download.pytorch.org/whl/cu118 torch torchvision torchaudio
```

# Install somatix in this conda environment:
pip install .                                    # installs in ~5 seconds
```

Alternatively, after activating the **somatix** environment, run the somatix runner directly:

```bash
python ./source/somatix-runner.py [OPTIONS]
```

SomatiX candidate extraction uses `allele_counter`, a C++ executable written
to speed up pileup-style allele counting from indexed BAM files. A precompiled
Linux binary is provided at `source/somatix/bin/allele_counter` for convenience,
and the source code is provided at `source/somatix/bin/allele_counter.cpp`.
Because the bundled binary is dynamically linked, it requires compatible system
and `htslib` shared libraries at runtime, especially `glibc`, `libstdc++` and
`libhts`. If the bundled binary fails to load because of missing or incompatible
shared libraries, rebuild `allele_counter` from source on the target system or
use the Singularity image. See
[docs/allele-counter-troubleshooting.md](docs/allele-counter-troubleshooting.md)
for build and library troubleshooting.

### Option 2: SingularityCE (recommended)
Ensure you have **SingularityCE** installed, then pull the `somatix` container
from the Sylabs Library to use it instantly without installation:

```bash
singularity pull library://zlliu95/somatix/somatix:latest
```

Use SingularityCE for this `library://` pull command. Apptainer is a separate
fork of Singularity and may use a different default remote; an Apptainer binary
exposed as `singularity` may not be able to pull from the Sylabs `library://`
endpoint.

The Singularity image embeds the required `allele_counter` binary and its
runtime library environment, so users do not need to rebuild `allele_counter` or
manage its shared-library compatibility separately. Also ensure your FASTA and
BAM files are indexed (`samtools faidx` and `samtools index`).

---

## Case Study

Complete case-study workflows, including example data download, SomatiX
calling, `som.py` benchmarking and candidate-restricted summary analysis, are
available for:

- ONT: [docs/somatix-case-study-ont.md](docs/somatix-case-study-ont.md)
- PacBio HiFi: [docs/somatix-case-study-pacbio.md](docs/somatix-case-study-pacbio.md)

---

## Quick Start and Overview

SomatiX provides a single CLI with these primary subcommands:

- `candidates` — extract candidate SNVs from a BAM (uses `allele_counter`)
- `features` — generate compact HDF5 shard feature files from a candidate table and BAM
- `predict` — predict somatic/germline classes from paired tumor/normal feature shards using a PyTorch model
- `perm` — experimental feature-testing command that predicts after permuting one selected feature group across samples
- `call` — run the full pipeline (candidates -> features -> predict)
- `clean` — remove intermediate files

The CLI entrypoint is the wrapper `source/somatix-runner.py` which calls `somatix/somatix.py`.

## Model Architecture

SomatiX uses paired tumor and matched-normal branches with the same feature
definitions and network architecture. In each branch, the sequence/base-fraction
tensor is processed by convolutional neural network layers, while coverage,
base-quality, mapping-quality, mismatch-rate and target-position coverage
features are incorporated through fully connected layers. The tumor and normal
branch outputs are concatenated and passed through fully connected layers to
produce four-class logits for reference, heterozygous germline, homozygous
germline and somatic candidates.

The current prediction architecture defines **382,720 trainable parameters**.
Of these, **382,612 parameters are used by the active forward path**; the
remaining 108 parameters belong to a defined raw target-coverage branch that is
retained in the module for checkpoint compatibility but is not concatenated into
the current classifier forward path.

## Pretrained Models

SomatiX uses platform-specific PyTorch checkpoints for prediction. The model
architecture is the same across checkpoints, but separate weights are provided
for each sequencing platform and training set. Place downloaded checkpoints
under the repository-level `model/` directory:

```text
SomatiX/
  model/
    somatix_ont_multicancer.pth
    somatix_ont_hcc1395.pth
    somatix_pacbio_multicancer.pth
    somatix_pacbio_hcc1395.pth
```

Available model types:

| Checkpoint | Platform | Training set | Recommended use |
| --- | --- | --- | --- |
| `somatix_ont_multicancer.pth` | ONT | Multi-cancer model trained across the cancer cell-line panel | Default ONT model for general tumor-normal ONT somatic SNV calling. |
| `somatix_ont_hcc1395.pth` | ONT | HCC1395-specific model | ONT HCC1395 analyses or controlled comparisons against the multi-cancer ONT model. |
| `somatix_pacbio_multicancer.pth` | PacBio HiFi | Multi-cancer model trained across the cancer cell-line panel | Default PacBio HiFi model for general tumor-normal PacBio somatic SNV calling. |
| `somatix_pacbio_hcc1395.pth` | PacBio HiFi | HCC1395-specific model | PacBio HCC1395 analyses or controlled comparisons against the multi-cancer PacBio model. |

Use the checkpoint that matches the sequencing platform of the input BAM files.
For most external samples, start with the corresponding multi-cancer model. The
HCC1395-specific checkpoints are mainly intended for HCC1395-only analyses and
benchmarking experiments that compare single-cell-line and multi-cancer model
training.

## Basic Usage

### 1) Extract Candidates

```bash
somatix candidates \
  --bam sample_case.bam \
  --ref reference.fa \
  --output case_candidates.txt \
  --allele-counter /path/to/allele_counter
```

Outputs a tab-delimited candidate table (columns include chrom, pos, ref, alt, ref_reads, alt_reads, total_coverage, vaf, ...).

### 2) Extract Tumor and Normal Features

```bash
# Case/tumor features at case candidate sites
somatix features \
  --candidates case_candidates.txt \
  --bam sample_case.bam \
  --ref reference.fa \
  --output-prefix case_features \
  --filtered-prefix case_candidates.filtered

# Normal features at the same tumor candidate sites
somatix features \
  --candidates case_candidates.txt \
  --bam sample_normal.bam \
  --ref reference.fa \
  --output-prefix control_features \
  --filtered-prefix case_candidates.filtered
```

The `features` subcommand processes one BAM per run. Run it once for the case/tumor BAM and once for the matched normal BAM, using the same `--candidates` and `--filtered-prefix` so both feature sets are aligned to the same variant loci. This generates directories like `case_features.chr1/shard_000000.h5`, `control_features.chr1/shard_000000.h5`, and `manifest.tsv` for each processed chromosome.

### 3) Predict from Features

```bash
somatix predict \
  --case-features-prefix case_features \
  --control-features-prefix control_features \
  --variant-prefix case_candidates.filtered \
  --model model/somatix_ont_multicancer.pth \
  --output somatix_predict.txt \
  --bam sample_case.bam
```

Produces `somatix_predict.txt`, split somatic/germline tables, and VCFs: `somatix_predict.vcf`, `somatix_predict.somatic.vcf`, and `somatix_predict.germline.vcf`.

### 4) Full Pipeline

```bash
somatix call \
  --bam-case sample_case.bam \
  --bam-control sample_normal.bam \
  --ref reference.fa \
  --model model/somatix_ont_multicancer.pth \
  --outdir ./somatix_out
```

- Use `--skip-prediction` to stop after generating features and skip the DL prediction step.

### 5) Clean Intermediate Files

```bash
somatix clean --outdir ./somatix_out
```

---

## Detailed CLI Options

The CLI entrypoint supports `--version` at the top level:

| Parameter | Required | Default | Description |
| --- | --- | --- | --- |
| `--version` | No | N/A | Print the installed SomatiX version and exit. |

### `somatix candidates`

Extract candidate SNVs from one BAM using the `allele_counter` backend.

| Parameter | Required | Default | Description |
| --- | --- | --- | --- |
| `--bam` | Yes | N/A | Input BAM file to scan for candidate SNVs. |
| `--ref` | Yes | N/A | Indexed reference FASTA file. |
| `--output` | Yes | N/A | Output tab-delimited candidate table. |
| `--region` | No | Whole genome | Optional genomic interval, for example `chr1:100000-200000`. |
| `--threads` | No | `8` | Number of threads for candidate extraction. |
| `--allele-counter` | No | `allele_counter` | Path to the `allele_counter` executable. |
| `--min-alt` | No | `3` | Minimum ALT-read count required for a candidate. |
| `--min-vaf`, `--vaf` | No | `0.05` | Minimum variant allele fraction required for a candidate. |
| `--min-total-coverage`, `--min-total` | No | `3` | Minimum total read coverage required for a candidate. |
| `--min-mapq`, `--min-MQ` | No | `20` | Minimum read mapping quality used during candidate extraction. |
| `--min-baseq`, `--min-BQ` | No | `10` | Minimum base quality used during candidate extraction. |
| `--max-depth` | No | `5000` | Maximum depth used during candidate extraction. |
| `--excl-flags` | No | `2316` | SAM flag bitmask for reads to exclude. |
| `--block-size` | No | `1000000` | Genomic block size used while scanning the BAM. |

### `somatix features`

Filter candidate variants and extract compact HDF5 feature shards from a BAM.

| Parameter | Required | Default | Description |
| --- | --- | --- | --- |
| `--candidates` | Yes | N/A | Input candidate table, usually from `somatix candidates`. |
| `--bam` | Yes | N/A | BAM file used to extract features at candidate loci. |
| `--ref` | Yes | N/A | Indexed reference FASTA file. |
| `--output-prefix` | Yes | N/A | Prefix for per-chromosome HDF5 shard directories, such as `case_features`. |
| `--filtered-prefix` | Yes | N/A | Prefix for filtered per-chromosome candidate tables. |
| `--threads` | No | `8` | Number of threads for feature extraction. |
| `--region` | No | Whole genome | Optional genomic interval, for example `chr1:100000-200000`. |
| `--feature-depth`, `--depth` | No | `1000` | Maximum reads sampled per candidate site during feature extraction. |
| `--chunk-bp` | No | `1000` | Genomic chunk size for chunk-based feature extraction. |
| `--min-vaf`, `--vaf` | No | `0.05` | Minimum VAF retained when filtering the candidate table. |
| `--min-total-coverage`, `--min-total` | No | `3` | Minimum total coverage retained when filtering the candidate table. |
| `--min-alt` | No | `3` | Minimum ALT-read count retained when filtering the candidate table. |
| `--shard-size` | No | `100000` | Maximum number of candidate records per HDF5 shard. |
| `--compression` | No | No compression | HDF5 compression method; allowed values are `lzf` and `gzip`. Omit this option to write uncompressed shards. |

### `somatix predict`

Predict somatic and germline classes from paired tumor/normal HDF5 feature shards.

| Parameter | Required | Default | Description |
| --- | --- | --- | --- |
| `--case-features-prefix` | Yes | N/A | Prefix for case feature shard directories, such as `case_features`. |
| `--control-features-prefix` | Yes | N/A | Prefix for matched-normal/control feature shard directories, such as `control_features`. |
| `--variant-prefix` | Yes | N/A | Prefix for filtered candidate tables used to align variants with feature shards. |
| `--model` | Yes | N/A | PyTorch model checkpoint used for prediction. |
| `--output` | Yes | N/A | Output prediction table path. |
| `--bam` | No | N/A | Optional BAM path passed to the prediction output/annotation step. |
| `--device` | No | `cpu` | Device for model inference, such as `cpu` or a CUDA device string. |
| `--region` | No | Whole genome | Optional genomic interval; prediction checks and runs only the selected chromosome. |

### `somatix perm`

Experimental feature-testing command, not the standard SomatiX prediction
workflow. It runs prediction with one feature group permuted across samples
before model inference. For grouped features, all related ref/alt and
forward/reverse arrays are permuted together using the same sample order. The
inputs and outputs match `somatix predict`: it writes the prediction table,
somatic/germline split tables, and VCF files using the same output naming
convention. Use the resulting VCFs for benchmark-based feature importance
analysis against a truth set.

| Parameter | Required | Default | Description |
| --- | --- | --- | --- |
| `--case-features-prefix` | Yes | N/A | Prefix for case feature shard directories, such as `case_features`. |
| `--control-features-prefix` | Yes | N/A | Prefix for matched-normal/control feature shard directories, such as `control_features`. |
| `--variant-prefix` | Yes | N/A | Prefix for filtered candidate tables used to align variants with feature shards. |
| `--model` | Yes | N/A | PyTorch model checkpoint used for prediction. |
| `--output` | Yes | N/A | Output prediction table path. VCFs are written beside this file using the same convention as `predict`. |
| `--feature` | Yes | N/A | One feature group to permute across samples. Allowed values are listed below. |
| `--sample` | No | `both` | Which sample tower to permute: `case`, `control`, or `both`. Here `control` refers to the matched normal sample tower. |
| `--bam` | No | N/A | Optional BAM path passed to the prediction output/annotation step for VCF contigs. |
| `--device` | No | `cpu` | Device for model inference, such as `cpu` or a CUDA device string. |
| `--seed` | No | `2026` | Random seed used to shuffle sample order. |
| `--region` | No | Whole genome | Optional genomic interval; prediction checks and runs only the selected chromosome. |

Allowed `--feature` group values:

```text
base_fraction
target_position_coverage
coverage
mapq
baseq
mismatch_rate
```

Feature group details:

| Feature group | Meaning | HDF5 datasets permuted together | Per-sample tensor structure |
| --- | --- | --- | --- |
| `base_fraction` | Five-channel base-fraction branch input: padded sequence context plus ref/alt A/C/G/T/DEL/INS base-fraction pileups on forward and reverse strands. | `one_hot_encoded_sequence`, `ref_base_fractions_forward`, `ref_base_fractions_reverse`, `alt_base_fractions_forward`, `alt_base_fractions_reverse` | One sequence tensor, shape `(4, 61)`, plus four base-fraction tensors, each shape `(6, 61)`. The model pads sequence to six rows internally, then stacks these five tensors as one branch. |
| `target_position_coverage` | Normalized total coverage at the candidate position. | `target_position_coverage` | One scalar tensor, shape `(1,)` |
| `coverage` | Ref/alt allele coverage pileups on forward and reverse strands. | `ref_coverage_forward`, `ref_coverage_reverse`, `alt_coverage_forward`, `alt_coverage_reverse` | Four tensors, each shape `(61,)` |
| `mapq` | Mean mapping quality for ref/alt-supporting reads on forward and reverse strands. | `ref_mapq_forward`, `ref_mapq_reverse`, `alt_mapq_forward`, `alt_mapq_reverse` | Four scalar tensors, each shape `(1,)` |
| `baseq` | Mean base quality for ref/alt-supporting reads on forward and reverse strands. | `ref_baseq_forward`, `ref_baseq_reverse`, `alt_baseq_forward`, `alt_baseq_reverse` | Four tensors, each shape `(61,)` |
| `mismatch_rate` | Mismatch-rate summaries for ref/alt-supporting reads on forward and reverse strands. | `ref_mismatch_rate_forward`, `ref_mismatch_rate_reverse`, `alt_mismatch_rate_forward`, `alt_mismatch_rate_reverse` | Four scalar tensors, each shape `(1,)` |

### `somatix call`

Run the full pipeline: case candidate extraction, paired tumor/normal feature extraction, and optional prediction.

| Parameter | Required | Default | Description |
| --- | --- | --- | --- |
| `--bam-case` | Yes | N/A | Case/tumor BAM file. |
| `--bam-control` | Yes | N/A | Matched normal BAM file. |
| `--ref` | Yes | N/A | Indexed reference FASTA file. |
| `--model` | Required unless `--skip-prediction` is set | N/A | PyTorch model checkpoint used for prediction. |
| `--outdir` | Yes | N/A | Output directory for candidates, features, and prediction results. |
| `--region` | No | Whole genome | Optional genomic interval, for example `chr1:100000-200000`. |
| `--threads` | No | `8` | Number of threads used by candidate and feature extraction. |
| `--allele-counter` | No | `allele_counter` | Path to the `allele_counter` executable. |
| `--device` | No | `cpu` | Device for model inference, such as `cpu` or a CUDA device string. |
| `--skip-prediction` | No | `False` | Stop after candidate and feature extraction; do not run DL prediction. |
| `--min-alt` | No | `3` | Minimum ALT-read count for candidate extraction and feature-input filtering. |
| `--min-vaf`, `--vaf` | No | `0.05` | Minimum VAF for candidate extraction and feature-input filtering. |
| `--min-total-coverage`, `--min-total` | No | `3` | Minimum total coverage for candidate extraction and feature-input filtering. |
| `--min-mapq`, `--min-MQ` | No | `20` | Minimum read mapping quality used during candidate extraction. |
| `--min-baseq`, `--min-BQ` | No | `10` | Minimum base quality used during candidate extraction. |
| `--max-depth` | No | `5000` | Maximum depth used during candidate extraction. |
| `--excl-flags` | No | `2316` | SAM flag bitmask for reads to exclude. |
| `--block-size` | No | `1000000` | Genomic block size used while scanning the BAM. |
| `--feature-depth`, `--depth` | No | `1000` | Maximum reads sampled per candidate site during feature extraction. |
| `--chunk-bp` | No | `1000` | Genomic chunk size for chunk-based feature extraction. |
| `--shard-size` | No | `100000` | Maximum number of candidate records per HDF5 shard. |
| `--compression` | No | No compression | HDF5 compression method; allowed values are `lzf` and `gzip`. Omit this option to write uncompressed shards. |

### `somatix clean`

Remove intermediate files generated by `somatix call`.

| Parameter | Required | Default | Description |
| --- | --- | --- | --- |
| `--outdir` | Yes | N/A | SomatiX output directory containing intermediate files to remove. |

## Outputs

- Tab-separated prediction table (annotated): `<output>.txt`
- Split tables: `<output>.somatic.txt`, `<output>.germline.txt`
- VCFs: `<output>.vcf`, `<output>.somatic.vcf`, `<output>.germline.vcf`
- HDF5 shards: `shard_*.h5` under `case_features.chr*` and `control_features.chr*`, plus `manifest.tsv`.

## Resource Notes

- Tuning `--min-alt`, `--depth`, and `--chunk-bp` affects memory/runtime. For whole-genome data, use many threads and adjust shard sizes to trade memory for I/O.

---

## License
This project is distributed under the **SomatiX Research Use License** for academic, research, and non-commercial use. Commercial use requires a separate written license from the copyright holder. See the [LICENSE](LICENSE) file for details.

---

## Contact
For questions, bug reports, or feature requests, contact:

**Zelin Liu** – zlliu95@outlook.com
