# Troubleshooting `allele_counter`

`allele_counter` is the C++ executable used by SomatiX candidate extraction. It
performs pileup-style allele counting from indexed BAM files and is called by
`somatix candidates` and `somatix call`.

SomatiX provides:

- a precompiled Linux binary: `source/somatix/bin/allele_counter`
- the source code: `source/somatix/bin/allele_counter.cpp`

If the precompiled binary does not run on your system, either use the
SingularityCE container or rebuild `allele_counter` locally against your own
system libraries.

## Recommended Option: Use SingularityCE

The SingularityCE image embeds `allele_counter` and the compatible runtime
libraries:

```bash
singularity pull library://zlliu95/somatix/somatix:latest
singularity exec somatix_latest.sif somatix --version
```

Use SingularityCE for the `library://` pull command. Apptainer is a separate
fork of Singularity and may use a different default remote.

## Check Whether the Bundled Binary Works

From the SomatiX repository root:

```bash
ALLELE_COUNTER="${PWD}/source/somatix/bin/allele_counter"

"${ALLELE_COUNTER}" --help
ldd "${ALLELE_COUNTER}"
```

If `ldd` reports `not found`, a required shared library is missing from your
runtime environment.

Common examples:

```text
libhts.so.3: cannot open shared object file: No such file or directory
GLIBC_2.38 not found
```

`GLIBC_* not found` means the binary was built on a newer Linux system than the
one used to run it. Rebuild locally or use the SingularityCE image.

## Install Build Dependencies

On Ubuntu/Debian systems, install a compiler and HTSlib development headers:

```bash
sudo apt update
sudo apt install -y \
  build-essential \
  g++ \
  zlib1g-dev \
  libbz2-dev \
  liblzma-dev \
  libcurl4-openssl-dev \
  libhts-dev
```

If `libhts-dev` is unavailable or too old, install HTSlib from source or a
Conda/Bioconda environment.

## Rebuild with System HTSlib

From the SomatiX repository root:

```bash
cd source/somatix/bin

g++ -O3 -std=c++17 allele_counter.cpp -o allele_counter \
  -lhts -lz -lbz2 -llzma -lcurl -lpthread

chmod +x allele_counter
./allele_counter --help
ldd ./allele_counter
```

Then use the rebuilt binary explicitly:

```bash
somatix candidates \
  --bam sample_case.bam \
  --ref reference.fa \
  --output case_candidates.txt \
  --allele-counter "${PWD}/source/somatix/bin/allele_counter"
```

## Rebuild with a Custom HTSlib Path

If HTSlib was installed outside the system library path, point the compiler to
the HTSlib include and library directories:

```bash
HTSLIB_DIR=/path/to/htslib

cd source/somatix/bin

g++ -O3 -std=c++17 allele_counter.cpp -o allele_counter \
  -I "${HTSLIB_DIR}" \
  -L "${HTSLIB_DIR}" \
  -lhts -lz -lbz2 -llzma -lcurl -lpthread
```

If the binary compiles but cannot find `libhts` at runtime, set
`LD_LIBRARY_PATH`:

```bash
export LD_LIBRARY_PATH="${HTSLIB_DIR}:${LD_LIBRARY_PATH:-}"
./allele_counter --help
```

You can also check exactly which library is being used:

```bash
ldd ./allele_counter | grep -E "hts|z|bz2|lzma|curl"
```

## Install the Rebuilt Binary

The simplest approach is to keep the rebuilt binary at:

```text
source/somatix/bin/allele_counter
```

Then pass it with `--allele-counter` when using a source checkout.

If you want it available globally:

```bash
sudo install -m 0755 source/somatix/bin/allele_counter /usr/local/bin/allele_counter
which allele_counter
allele_counter --help
```

## Common Errors

### `allele_counter: command not found`

The binary is not in `PATH`. Use an explicit path:

```bash
--allele-counter /absolute/path/to/allele_counter
```

or install it into a directory in `PATH`, such as `/usr/local/bin`.

### `libhts.so.*: cannot open shared object file`

HTSlib is missing or the runtime linker cannot find it.

Options:

- install `libhts-dev` / `libhts`
- rebuild `allele_counter` against the HTSlib on your system
- set `LD_LIBRARY_PATH` to the directory containing `libhts.so`
- use the SingularityCE image

### HTSlib is installed, but the requested `libhts.so.*` file does not exist

Sometimes HTSlib is installed successfully, but `allele_counter` still fails
with an error like:

```text
libhts.so.3: cannot open shared object file: No such file or directory
```

This can happen when the bundled `allele_counter` binary was linked against one
versioned HTSlib shared-library name, but the system provides a different one.
For example, the binary may request `libhts.so.3`, while your installed HTSlib
directory contains files such as `libhts.so`, `libhts.so.1.22` or
`libhts.so.1.23.1`.

Check what the binary requests:

```bash
ldd source/somatix/bin/allele_counter | grep hts
```

Check which HTSlib files are available:

```bash
find /usr /usr/local "$CONDA_PREFIX" -name 'libhts.so*' 2>/dev/null
```

The safest fix is to rebuild `allele_counter` on the same system where it will
be run, so the binary links to the HTSlib version actually installed there.

If you know the installed HTSlib is ABI-compatible, a local symlink can also
solve the loader error. For example:

```bash
mkdir -p "${PWD}/local_lib"
ln -s /path/to/libhts.so.1.23.1 "${PWD}/local_lib/libhts.so.3"
export LD_LIBRARY_PATH="${PWD}/local_lib:${LD_LIBRARY_PATH:-}"
ldd source/somatix/bin/allele_counter | grep hts
```

Use this symlink approach carefully. It only works when the installed HTSlib is
compatible with the binary. If the binary starts but behaves unexpectedly,
remove the symlink and rebuild `allele_counter` from source instead.

### `GLIBC_2.xx not found`

The bundled binary was built on a newer Linux system. Rebuild locally on the
target machine or use the SingularityCE image.

### `undefined reference to hts_*` during compilation

The linker cannot find HTSlib, or the `-lhts` argument is missing or ordered
incorrectly. Keep source files before library flags:

```bash
g++ -O3 -std=c++17 allele_counter.cpp -o allele_counter \
  -lhts -lz -lbz2 -llzma -lcurl -lpthread
```

If HTSlib is in a custom location, add `-I` and `-L` paths as shown above.

## Quick Functional Test

Run on a small genomic interval:

```bash
allele_counter \
  --bam sample_case.bam \
  --ref reference.fa \
  --region chr1:1-1000000 \
  --min-mapq 20 \
  --min-baseq 10 \
  --min-alt 3 \
  --min-total-coverage 3 \
  --min-vaf 0.05 \
  --max-depth 5000
```

The BAM and FASTA must be indexed:

```bash
samtools index sample_case.bam
samtools faidx reference.fa
```
