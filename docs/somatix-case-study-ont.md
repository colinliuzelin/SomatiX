# SomatiX ONT Case Study

This case study shows an example SomatiX workflow for matched tumor-control
Oxford Nanopore Technologies (ONT) data. The example uses HCC1395 chromosome 1
so that the commands can be run quickly and compared against the SEQC2 somatic
truth set.

This document is for ONT. A separate PacBio HiFi case-study document can be
added later using the same structure.

The commands below assume that you run the example from the SomatiX repository
root:

```bash
cd /mnt/windows/mydata/data1/DNA_somatic/git/test_shard_chunk/SomatiX
```

## Data Details

This case study uses:

- matched HCC1395 tumor and normal ONT BAM files restricted to chromosome 1;
- a GRCh38 chromosome 1 reference FASTA;
- the SEQC2 HCC1395 high-confidence BED and somatic SNV truth VCF;
- a pretrained SomatiX ONT model checkpoint supplied by the user.

SomatiX evaluates candidate SNVs from the tumor BAM and extracts paired tumor
and control features at the same candidate sites. By default, the example uses
the candidate definition used in the SomatiX manuscript analyses:

- ALT reads at least 3;
- VAF at least 0.05;
- read mapping quality at least 20;
- base quality at least 10;
- SNV candidates only.

## Prepare Environment

### Option 1: Conda

```bash
git clone <SomatiX GitHub URL> SomatiX
cd SomatiX

conda env create -f environment.yml -n somatix
conda activate somatix
pip install .
```

The examples below call the runner script directly:

```bash
SOMATIX_DIR="${PWD}"
PYTHON_BIN="$(which python)"
```

If SomatiX is installed as a command-line entry point in your environment, you
can replace:

```bash
${PYTHON_BIN} ${SOMATIX_DIR}/source/somatix-runner.py
```

with:

```bash
somatix
```

### Option 2: Singularity

```bash
singularity pull library://zelinliu/somatix/somatix:latest
SOMATIX_SIF="${PWD}/somatix_latest.sif"
```

The example was tested with Apptainer/Singularity 1.1.9. If a local
installation has problems converting Docker images to SIF format, the following
commands install Apptainer 1.1.9 on Ubuntu/Debian systems:

```bash
cd /tmp
wget https://github.com/apptainer/apptainer/releases/download/v1.1.9/apptainer_1.1.9_amd64.deb
sudo apt remove -y singularity apptainer
sudo dpkg -i apptainer_1.1.9_amd64.deb
sudo apt -f install -y
```

The commands below are written for a Conda/source checkout. To run the same
commands with Singularity, prepend the command with:

```bash
singularity exec \
  -B "${BASE}:${BASE}" \
  -B "$(dirname "${MODEL}")":"$(dirname "${MODEL}")" \
  -B "$(dirname "${ALLELE_COUNTER}")":"$(dirname "${ALLELE_COUNTER}")" \
  "${SOMATIX_SIF}"
```

## Download Example Input Data

```bash
BASE="${PWD}/example/ont"
INPUT_DIR="${BASE}/input"
OUTPUT_DIR="${BASE}/output"
MODEL_DIR="${SOMATIX_DIR}/model"

mkdir -p "${INPUT_DIR}" "${OUTPUT_DIR}"

HTTPDIR="https://storage.googleapis.com/deepvariant/deepsomatic-case-studies/deepsomatic-chr1-case-studies"
SEQC2_HTTPDIR="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/release/latest"

# GRCh38 chromosome 1 reference.
wget -O "${INPUT_DIR}/GRCh38.chr1.fa" \
  "${HTTPDIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.chr1.fna"
wget -O "${INPUT_DIR}/GRCh38.chr1.fa.fai" \
  "${HTTPDIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.chr1.fna.fai"

# HCC1395 ONT tumor-control BAM files.
wget -O "${INPUT_DIR}/HCC1395_ont.normal.chr1.bam" \
  "${HTTPDIR}/HCC1395_ont.normal.chr1.bam"
wget -O "${INPUT_DIR}/HCC1395_ont.normal.chr1.bam.bai" \
  "${HTTPDIR}/HCC1395_ont.normal.chr1.bam.bai"
wget -O "${INPUT_DIR}/HCC1395_ont.tumor.chr1.bam" \
  "${HTTPDIR}/HCC1395_ont.tumor.chr1.bam"
wget -O "${INPUT_DIR}/HCC1395_ont.tumor.chr1.bam.bai" \
  "${HTTPDIR}/HCC1395_ont.tumor.chr1.bam.bai"

# SEQC2 HCC1395 benchmark resources.
wget -O "${INPUT_DIR}/high-confidence_sSNV_in_HC_regions_v1.2.vcf.gz" \
  "${SEQC2_HTTPDIR}/high-confidence_sSNV_in_HC_regions_v1.2.vcf.gz"
wget -O "${INPUT_DIR}/High-Confidence_Regions_v1.2.bed" \
  "${SEQC2_HTTPDIR}/High-Confidence_Regions_v1.2.bed"
```

Set paths used by the remaining commands:

```bash
REF="${INPUT_DIR}/GRCh38.chr1.fa"
TUMOR_BAM="${INPUT_DIR}/HCC1395_ont.tumor.chr1.bam"
NORMAL_BAM="${INPUT_DIR}/HCC1395_ont.normal.chr1.bam"
TRUTH_VCF="${INPUT_DIR}/high-confidence_sSNV_in_HC_regions_v1.2.vcf.gz"
TRUTH_BED="${INPUT_DIR}/High-Confidence_Regions_v1.2.bed"

# Provide one or both pretrained SomatiX ONT model checkpoints.
# The multi-cancer checkpoint is the default model for this case study.
MODEL_MULTICANCER="${MODEL_DIR}/somatix_ont_multicancer.pth"
MODEL_HCC1395="${MODEL_DIR}/somatix_ont_hcc1395.pth"
MODEL="${MODEL_MULTICANCER}"

# Use the bundled allele_counter binary by default.
ALLELE_COUNTER="${SOMATIX_DIR}/source/somatix/bin/allele_counter"
```

## Run SomatiX With One Command

The `call` subcommand runs candidate extraction, paired tumor/control feature
extraction and prediction.

```bash
${PYTHON_BIN} "${SOMATIX_DIR}/source/somatix-runner.py" call \
  --bam-case "${TUMOR_BAM}" \
  --bam-control "${NORMAL_BAM}" \
  --ref "${REF}" \
  --model "${MODEL}" \
  --allele-counter "${ALLELE_COUNTER}" \
  --outdir "${OUTPUT_DIR}/somatix_chr1" \
  --region chr1 \
  --threads "$(nproc)" \
  --min-baseq 10 \
  --min-mapq 20 \
  --vaf 0.05 \
  --min-alt 3 \
  --min-total 3 \
  --max-depth 5000 \
  --shard-size 100000 \
  --chunk-bp 1000
```

Expected main outputs:

```text
${OUTPUT_DIR}/somatix_chr1/case_candidates.chr1
${OUTPUT_DIR}/somatix_chr1/case_candidates.filtered.chr1
${OUTPUT_DIR}/somatix_chr1/case_features.chr1/
${OUTPUT_DIR}/somatix_chr1/control_features.chr1/
${OUTPUT_DIR}/somatix_chr1/somatix_predict.txt
${OUTPUT_DIR}/somatix_chr1/somatix_predict.vcf
${OUTPUT_DIR}/somatix_chr1/somatix_predict.somatic.vcf
${OUTPUT_DIR}/somatix_chr1/somatix_predict.germline.vcf
```

The standard VCF used for somatic benchmarking is:

```bash
QUERY_VCF="${OUTPUT_DIR}/somatix_chr1/somatix_predict.vcf"
```

## Run SomatiX Step By Step

The full pipeline can also be run as separate candidate, feature and prediction
steps. This is useful for debugging, reusing feature shards or running
permutation feature tests.

### 1. Candidate SNV Extraction

```bash
${PYTHON_BIN} "${SOMATIX_DIR}/source/somatix-runner.py" candidates \
  --bam "${TUMOR_BAM}" \
  --ref "${REF}" \
  --output "${OUTPUT_DIR}/case_candidates.txt" \
  --region chr1 \
  --threads "$(nproc)" \
  --allele-counter "${ALLELE_COUNTER}" \
  --min-baseq 10 \
  --min-mapq 20 \
  --vaf 0.05 \
  --min-alt 3 \
  --min-total 3 \
  --max-depth 5000
```

### 2. Tumor Feature Extraction

```bash
${PYTHON_BIN} "${SOMATIX_DIR}/source/somatix-runner.py" features \
  --candidates "${OUTPUT_DIR}/case_candidates.txt" \
  --bam "${TUMOR_BAM}" \
  --ref "${REF}" \
  --output-prefix "${OUTPUT_DIR}/case_features" \
  --filtered-prefix "${OUTPUT_DIR}/case_candidates.filtered" \
  --region chr1 \
  --threads "$(nproc)" \
  --min-baseq 10 \
  --min-mapq 20 \
  --vaf 0.05 \
  --min-alt 3 \
  --min-total 3 \
  --shard-size 100000 \
  --chunk-bp 1000
```

### 3. Control Feature Extraction

The control BAM is processed at the same tumor candidate loci. Use the same
`--candidates` and `--filtered-prefix` paths so tumor and control features stay
aligned.

```bash
${PYTHON_BIN} "${SOMATIX_DIR}/source/somatix-runner.py" features \
  --candidates "${OUTPUT_DIR}/case_candidates.txt" \
  --bam "${NORMAL_BAM}" \
  --ref "${REF}" \
  --output-prefix "${OUTPUT_DIR}/control_features" \
  --filtered-prefix "${OUTPUT_DIR}/case_candidates.filtered" \
  --region chr1 \
  --threads "$(nproc)" \
  --min-baseq 10 \
  --min-mapq 20 \
  --vaf 0.05 \
  --min-alt 3 \
  --min-total 3 \
  --shard-size 100000 \
  --chunk-bp 1000
```

### 4. Prediction

```bash
${PYTHON_BIN} "${SOMATIX_DIR}/source/somatix-runner.py" predict \
  --case-features-prefix "${OUTPUT_DIR}/case_features" \
  --control-features-prefix "${OUTPUT_DIR}/control_features" \
  --variant-prefix "${OUTPUT_DIR}/case_candidates.filtered" \
  --model "${MODEL}" \
  --output "${OUTPUT_DIR}/somatix_predict.txt" \
  --bam "${TUMOR_BAM}" \
  --region chr1
```

The prediction VCF is written next to the prediction table:

```bash
QUERY_VCF="${OUTPUT_DIR}/somatix_predict.vcf"
```

## Benchmark With som.py

The following example benchmarks the SomatiX VCF against the SEQC2 HCC1395
truth set using `som.py` from hap.py. Only chromosome 1 is evaluated.

```bash
HAP_PY_SIF="${PWD}/hap.py_latest.sif"
singularity pull "${HAP_PY_SIF}" docker://pkrusche/hap.py:latest

SOMPY_OUT="${OUTPUT_DIR}/sompy_output/somatix_chr1"
mkdir -p "$(dirname "${SOMPY_OUT}")"

singularity exec \
  -B "${INPUT_DIR}:${INPUT_DIR}" \
  -B "${OUTPUT_DIR}:${OUTPUT_DIR}" \
  "${HAP_PY_SIF}" \
  /opt/hap.py/bin/som.py \
  -N "${TRUTH_VCF}" \
  "${QUERY_VCF}" \
  -r "${REF}" \
  -o "${SOMPY_OUT}" \
  --feature-table generic \
  -R "${TRUTH_BED}" \
  -l chr1
```

The main summary table is written to:

```text
${SOMPY_OUT}.stats.csv
```

The feature table used for downstream candidate-restricted analyses is written
to:

```text
${SOMPY_OUT}.features.csv
```

## Candidate-Restricted Summary

For analyses that match the SomatiX candidate definition, restrict the `som.py`
feature table to sites present in the filtered candidate table, using
`chrom:pos` as the site key. In the manuscript analyses, downstream SNV metrics
were restricted to sites with at least three ALT reads and VAF at least 0.05
from reads with mapping quality at least 20 and base quality at least 10.

Use the helper script below to calculate TP, FP, FN, precision, recall and F1
from the filtered `som.py` feature-table tags:

```bash
python "${SOMATIX_DIR}/other_scripts/benchmark/candidate_restricted_sompy_metrics.py" \
  --sompy-features "${SOMPY_OUT}.features.csv" \
  --candidates "${OUTPUT_DIR}/somatix_chr1/case_candidates.filtered.chr1" \
  --output "${OUTPUT_DIR}/sompy_output/somatix_chr1.candidate_restricted_metrics.csv" \
  --filtered-features-output "${OUTPUT_DIR}/sompy_output/somatix_chr1.candidate_restricted_features.csv" \
  --sample-id HCC1395 \
  --platform ONT \
  --tool SomatiX \
  --data-type raw \
  --chrom chr1 \
  --min-alt 3 \
  --min-vaf 0.05
```

## Clean Intermediate Files

If you only need the final prediction outputs, remove intermediate candidate
and feature files with:

```bash
${PYTHON_BIN} "${SOMATIX_DIR}/source/somatix-runner.py" clean \
  --outdir "${OUTPUT_DIR}/somatix_chr1"
```

## Notes

- Ensure the reference FASTA has a `.fai` index and each BAM has a `.bai` index.
- `--region chr1` is used here for a small case study. Remove this option or
  replace it with another interval for larger analyses.
- `--chunk-bp 1000` is the current default and is shown explicitly for
  reproducibility.
- Use `--device cuda:0` in `predict` or `call` if a compatible GPU and PyTorch
  installation are available.
