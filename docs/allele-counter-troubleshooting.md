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

## Install HTSlib and Build Dependencies

`allele_counter` directly depends on HTSlib. When linking against HTSlib, the
final executable may also need HTSlib's compression/network/threading
dependencies:

- `zlib`: linked with `-lz`
- `bzip2`: linked with `-lbz2`
- `xz/liblzma`: linked with `-llzma`
- `libcurl`: linked with `-lcurl`
- pthreads: linked with `-lpthread`

### Option A: Use the System HTSlib Package

On Ubuntu/Debian systems, install HTSlib and the required compiler/development
libraries with:

```bash
sudo apt update
sudo apt install -y \
  build-essential \
  g++ \
  libhts-dev \
  zlib1g-dev \
  libbz2-dev \
  liblzma-dev \
  libcurl4-openssl-dev
```

### Option B: Use Conda/Bioconda HTSlib

With Conda/Bioconda, install HTSlib and the related libraries with:

```bash
conda install -y -c conda-forge -c bioconda \
  htslib \
  zlib \
  bzip2 \
  xz \
  libcurl \
  cxx-compiler \
  make
```

If using the Conda/Bioconda HTSlib, compile `allele_counter` against the
headers and libraries in the active conda environment:

```bash
cd source/somatix/bin

g++ -O3 -std=c++17 allele_counter.cpp -o allele_counter \
  -I"${CONDA_PREFIX}/include" \
  -L"${CONDA_PREFIX}/lib" \
  -Wl,-rpath,"${CONDA_PREFIX}/lib" \
  -lhts -lz -lbz2 -llzma -lcurl -lpthread
```

### Option C: Build HTSlib From Source

If the system `libhts-dev` package is unavailable or too old, build HTSlib from
the official [samtools/htslib](https://github.com/samtools/htslib) source
repository. First install the compiler and libraries needed to build HTSlib:

```bash
sudo apt update
sudo apt install -y \
  build-essential \
  autoconf \
  automake \
  libtool \
  g++ \
  zlib1g-dev \
  libbz2-dev \
  liblzma-dev \
  libcurl4-openssl-dev
```

Then build HTSlib:

```bash
git clone https://github.com/samtools/htslib.git
cd htslib
autoreconf -i
./configure --prefix="${HOME}/.local"
make
make install
```

## Makefile Options

The `source/somatix/bin/Makefile` supports these variables:

| Variable | Default | Meaning |
| --- | --- | --- |
| `HTSLIB_DIR` | empty | Optional directory containing HTSlib headers and libraries. |
| `USE_RPATH` | `0` | Set to `1` to embed `HTSLIB_DIR` as a runtime library search path. |
| `PREFIX` | `/usr/local` | Installation prefix used by `make install`. |
| `CXX` | `g++` | C++ compiler. |
| `CXXFLAGS` | `-O3 -std=c++17` | C++ compiler flags. |

## Rebuild with System HTSlib

From the SomatiX repository root:

```bash
cd source/somatix/bin

make clean
make
./allele_counter --help
ldd ./allele_counter
```

Then use the rebuilt binary explicitly:

```bash
somatix candidates \
  --bam sample_case.bam \
  --ref reference.fa \
  --output case_candidates.txt \
  --allele-counter "${PWD}/allele_counter"
```

## Rebuild with a Custom HTSlib Path

If HTSlib was installed outside the system library path, point the compiler to
the HTSlib include and library directories:

```bash
HTSLIB_DIR=/path/to/htslib

cd source/somatix/bin

make clean
make HTSLIB_DIR="${HTSLIB_DIR}"
```

If the binary compiles but cannot find `libhts` at runtime, set
`LD_LIBRARY_PATH`:

```bash
export LD_LIBRARY_PATH="${HTSLIB_DIR}:${LD_LIBRARY_PATH:-}"
./allele_counter --help
```

Alternatively, build with an embedded runtime search path so the binary uses the
same HTSlib directory used during compilation:

```bash
make clean
make HTSLIB_DIR="${HTSLIB_DIR}" USE_RPATH=1
ldd ./allele_counter | grep hts
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
cd source/somatix/bin
sudo make install PREFIX=/usr/local
which allele_counter
allele_counter --help
```

To install somewhere else:

```bash
make install PREFIX="${HOME}/.local"
export PATH="${HOME}/.local/bin:${PATH}"
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

### Other shared libraries are missing

`allele_counter` may also fail if an HTSlib dependency is missing from the
runtime environment, for example:

```text
error while loading shared libraries: libbz2.so.1.0: cannot open shared object file: No such file or directory
```

Similar errors can occur for `libz`, `liblzma`, `libcurl`, `libstdc++` or other
shared libraries. Check missing libraries with:

```bash
ldd source/somatix/bin/allele_counter
```

On Ubuntu/Debian systems, install the common runtime/development packages with:

```bash
sudo apt update
sudo apt install -y \
  zlib1g \
  libbz2-1.0 \
  liblzma5 \
  libcurl4 \
  libstdc++6
```

With Conda/Bioconda, install the corresponding libraries in the active
environment and make sure the conda library directory is visible:

```bash
conda install -y -c conda-forge \
  zlib \
  bzip2 \
  xz \
  libcurl \
  libstdcxx-ng

export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"
ldd source/somatix/bin/allele_counter
```

If several libraries are missing, rebuilding `allele_counter` on the target
system or using the SingularityCE image is usually safer than manually fixing
one library at a time.

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
For Conda/Bioconda HTSlib, use the `${CONDA_PREFIX}` compile command shown in
the Conda section.

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
