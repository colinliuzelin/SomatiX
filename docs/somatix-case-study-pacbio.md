# SomatiX PacBio HiFi Case Study

This case study shows an example SomatiX workflow for matched tumor-normal
Pacific Biosciences (PacBio) HiFi data. The example uses HCC1395 chromosome 1
so that the commands can be run quickly and compared against the SEQC2 somatic
truth set.

For convenience and reproducibility, this example uses the same HCC1395 PacBio
HiFi case-study BAM files made available through the DeepSomatic GitHub
case-study workflow. We thank the DeepSomatic team for providing rich public
resources that support somatic variant-calling benchmarking.

This document is for PacBio HiFi. The ONT version is available at
[somatix-case-study-ont.md](somatix-case-study-ont.md).

The commands below assume that you run the example from the root path of the
SomatiX repository.

## Data Details

This case study uses:

- matched HCC1395 tumor and normal PacBio HiFi BAM files restricted to chromosome 1;
- a GRCh38 chromosome 1 reference FASTA;
- the SEQC2 HCC1395 high-confidence BED and somatic SNV truth VCF;
- a pretrained SomatiX PacBio HiFi model checkpoint supplied by the user.

SomatiX evaluates candidate SNVs from the tumor BAM and extracts paired tumor
and normal features at the same candidate sites. By default, the example uses
the candidate definition used in the SomatiX manuscript analyses:

- ALT reads at least 3;
- VAF at least 0.05;
- read mapping quality at least 20;
- base quality at least 10;
- SNV candidates only.

## Prepare Environment

### Option 1: Conda

If you have not cloned SomatiX yet, clone the repository first:

```bash
git clone <SomatiX GitHub URL> SomatiX
cd SomatiX
```

Then create and activate the Conda environment:

```bash
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

### Option 2: SingularityCE

```bash
singularity pull library://zlliu95/somatix/somatix:latest
SOMATIX_SIF="${PWD}/somatix_latest.sif"
```

Use SingularityCE for the `library://` pull command above. Apptainer is a
separate fork of Singularity and may use a different default remote; an
Apptainer binary exposed as `singularity` may not be able to pull from the
Sylabs `library://` endpoint.

The Singularity image includes the required `allele_counter` binary. The
container command below uses the bundled default and does not require an
external `--allele-counter` path.

## Download Example Input Data

```bash
BASE="${PWD}/example/pacbio"
INPUT_DIR="${BASE}/input"
OUTPUT_DIR="${BASE}/output"
SOMATIX_DIR="${PWD}"
MODEL_DIR="${SOMATIX_DIR}/model"

mkdir -p "${INPUT_DIR}" "${OUTPUT_DIR}"

HTTPDIR="https://storage.googleapis.com/deepvariant/deepsomatic-case-studies/deepsomatic-chr1-case-studies"
SEQC2_HTTPDIR="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/release/latest"

# GRCh38 chromosome 1 reference.
wget -O "${INPUT_DIR}/GRCh38.chr1.fa" \
  "${HTTPDIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.chr1.fna"
wget -O "${INPUT_DIR}/GRCh38.chr1.fa.fai" \
  "${HTTPDIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.chr1.fna.fai"

# HCC1395 PacBio HiFi tumor-normal BAM files.
wget -O "${INPUT_DIR}/HCC1395_pacbio.normal.chr1.bam" \
  "${HTTPDIR}/HCC1395_pacbio.normal.chr1.bam"
wget -O "${INPUT_DIR}/HCC1395_pacbio.normal.chr1.bam.bai" \
  "${HTTPDIR}/HCC1395_pacbio.normal.chr1.bam.bai"
wget -O "${INPUT_DIR}/HCC1395_pacbio.tumor.chr1.bam" \
  "${HTTPDIR}/HCC1395_pacbio.tumor.chr1.bam"
wget -O "${INPUT_DIR}/HCC1395_pacbio.tumor.chr1.bam.bai" \
  "${HTTPDIR}/HCC1395_pacbio.tumor.chr1.bam.bai"

# SEQC2 HCC1395 benchmark resources.
wget -O "${INPUT_DIR}/high-confidence_sSNV_in_HC_regions_v1.2.vcf.gz" \
  "${SEQC2_HTTPDIR}/high-confidence_sSNV_in_HC_regions_v1.2.vcf.gz"
wget -O "${INPUT_DIR}/High-Confidence_Regions_v1.2.bed" \
  "${SEQC2_HTTPDIR}/High-Confidence_Regions_v1.2.bed"
```

Set paths used by the remaining commands:

```bash
REF="${INPUT_DIR}/GRCh38.chr1.fa"
TUMOR_BAM="${INPUT_DIR}/HCC1395_pacbio.tumor.chr1.bam"
NORMAL_BAM="${INPUT_DIR}/HCC1395_pacbio.normal.chr1.bam"
TRUTH_VCF="${INPUT_DIR}/high-confidence_sSNV_in_HC_regions_v1.2.vcf.gz"
TRUTH_BED="${INPUT_DIR}/High-Confidence_Regions_v1.2.bed"

# Provide one or both pretrained SomatiX PacBio HiFi model checkpoints.
# The multi-cancer checkpoint is the default model for this case study.
MODEL_MULTICANCER="${MODEL_DIR}/somatix_pacbio_multicancer.pth"
MODEL_HCC1395="${MODEL_DIR}/somatix_pacbio_hcc1395.pth"
MODEL="${MODEL_MULTICANCER}"

# For a source/Conda run, use the allele_counter binary included in this checkout.
ALLELE_COUNTER="${SOMATIX_DIR}/source/somatix/bin/allele_counter"
```

## Run SomatiX with Conda

The `call` subcommand runs candidate extraction, paired tumor/normal feature
extraction and prediction. In a source/Conda checkout, use the local runner and
the local `allele_counter` binary:

```bash
${PYTHON_BIN} "${SOMATIX_DIR}/source/somatix-runner.py" call \
  --bam-case "${TUMOR_BAM}" \
  --bam-control "${NORMAL_BAM}" \
  --ref "${REF}" \
  --model "${MODEL}" \
  --allele-counter "${ALLELE_COUNTER}" \
  --outdir "${OUTPUT_DIR}/somatix_pacbio_chr1" \
  --region chr1 \
  --threads 24 \
  --min-baseq 10 \
  --min-mapq 20 \
  --vaf 0.05 \
  --min-alt 3 \
  --min-total 3 \
  --max-depth 5000 \
  --shard-size 100000 \
  --chunk-bp 1000
```

## Run SomatiX with Singularity

When using the Singularity image, run the same `call` workflow inside the
container. The image includes `allele_counter`, so no external
`--allele-counter` argument is needed:

```bash
singularity exec \
  -B "${BASE}:${BASE}" \
  -B "${MODEL_DIR}:${MODEL_DIR}" \
  "${SOMATIX_SIF}" \
  somatix call \
  --bam-case "${TUMOR_BAM}" \
  --bam-control "${NORMAL_BAM}" \
  --ref "${REF}" \
  --model "${MODEL}" \
  --outdir "${OUTPUT_DIR}/somatix_pacbio_chr1" \
  --region chr1 \
  --threads 24 \
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
${OUTPUT_DIR}/somatix_pacbio_chr1/case_candidates.chr1
${OUTPUT_DIR}/somatix_pacbio_chr1/case_candidates.filtered.chr1
${OUTPUT_DIR}/somatix_pacbio_chr1/case_features.chr1/
${OUTPUT_DIR}/somatix_pacbio_chr1/control_features.chr1/
${OUTPUT_DIR}/somatix_pacbio_chr1/somatix_predict.txt
${OUTPUT_DIR}/somatix_pacbio_chr1/somatix_predict.vcf
${OUTPUT_DIR}/somatix_pacbio_chr1/somatix_predict.somatic.vcf
${OUTPUT_DIR}/somatix_pacbio_chr1/somatix_predict.germline.vcf
```

The standard VCF used for somatic benchmarking is:

```bash
QUERY_VCF="${OUTPUT_DIR}/somatix_pacbio_chr1/somatix_predict.vcf"
```

## Benchmark with som.py

The following example benchmarks the SomatiX VCF against the SEQC2 HCC1395
truth set using `som.py` from hap.py. Only chromosome 1 is evaluated.

```bash
HAP_PY_SIF="${PWD}/hap.py_latest.sif"
singularity pull "${HAP_PY_SIF}" library://zlliu95/hap.py/hap.py:latest

SOMPY_OUT="${OUTPUT_DIR}/sompy_output/somatix_pacbio_chr1"
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

An example PacBio `som.py` summary table is shown below:

<table>
  <thead>
    <tr>
      <th>type</th><th>total.truth</th><th>total.query</th><th>tp</th><th>fp</th><th>fn</th><th>unk</th><th>ambi</th><th>recall</th><th>recall_lower</th><th>recall_upper</th><th>recall2</th><th>precision</th><th>precision_lower</th><th>precision_upper</th><th>na</th><th>ambiguous</th><th>fp.region.size</th><th>fp.rate</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>SNVs</td><td>3440</td><td>3337</td><td>3203</td><td>134</td><td>237</td><td>0</td><td>0</td><td>0.931105</td><td>0.922277</td><td>0.939208</td><td>0.931105</td><td>0.959844</td><td>0.95278</td><td>0.966111</td><td>0</td><td>0</td><td>248956422</td><td>0.538247</td>
    </tr>
    <tr>
      <td>records</td><td>3440</td><td>3337</td><td>3203</td><td>134</td><td>237</td><td>0</td><td>0</td><td>0.931105</td><td>0.922277</td><td>0.939208</td><td>0.931105</td><td>0.959844</td><td>0.95278</td><td>0.966111</td><td>0</td><td>0</td><td>248956422</td><td>0.538247</td>
    </tr>
  </tbody>
</table>

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
from the filtered `som.py` feature-table tags. Run this command in the Conda
environment used for the source workflow, or any Python environment that has the
script requirements installed:

```bash
python "${SOMATIX_DIR}/other_scripts/benchmark/candidate_restricted_sompy_metrics.py" \
  --sompy-features "${SOMPY_OUT}.features.csv" \
  --candidates "${OUTPUT_DIR}/somatix_pacbio_chr1/case_candidates.filtered.chr1" \
  --output "${OUTPUT_DIR}/sompy_output/somatix_pacbio_chr1.candidate_restricted_metrics.csv" \
  --filtered-features-output "${OUTPUT_DIR}/sompy_output/somatix_pacbio_chr1.candidate_restricted_features.csv" \
  --sample-id HCC1395 \
  --platform PacBio \
  --tool SomatiX \
  --data-type raw \
  --chrom chr1 \
  --min-alt 3 \
  --min-vaf 0.05
```

An example PacBio candidate-restricted summary table is shown below:

| sample_id | platform | tool | chrom | min_alt | min_vaf | candidate_sites | som.py SNV rows | rows in candidates | total.truth | total.query | tp | fp | fn | precision | recall | f1 |
| --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| HCC1395 | PacBio | SomatiX | chr1 | 3 | 0.05 | 316487 | 3574 | 3340 | 3206 | 3337 | 3203 | 134 | 3 | 0.959844 | 0.999064 | 0.979062 |

## Clean Intermediate Files

If you only need the final prediction outputs, remove intermediate candidate
and feature files.

For Conda/source runs:

```bash
${PYTHON_BIN} "${SOMATIX_DIR}/source/somatix-runner.py" clean \
  --outdir "${OUTPUT_DIR}/somatix_pacbio_chr1"
```

For Singularity runs:

```bash
singularity exec \
  -B "${BASE}:${BASE}" \
  "${SOMATIX_SIF}" \
  somatix clean \
  --outdir "${OUTPUT_DIR}/somatix_pacbio_chr1"
```

## Notes

- Ensure the reference FASTA has a `.fai` index and each BAM has a `.bai` index.
- `--region chr1` is used here for a small case study. Remove this option or
  replace it with another interval for larger analyses.
- `--chunk-bp 1000` is the current default and is shown explicitly for
  reproducibility.
- Use `--device cuda:0` in `predict` or `call` if a compatible GPU and PyTorch
  installation are available.
