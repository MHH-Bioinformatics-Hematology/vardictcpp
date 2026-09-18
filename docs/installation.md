# Installation

vardictcpp runs on Linux and macOS, on both x86-64 and ARM (Apple Silicon). Its only runtime
dependency is [htslib](https://github.com/samtools/htslib).

## Bioconda (recommended)

```bash
conda install -c bioconda vardictcpp
```

This pulls in a matching htslib automatically. [Mamba](https://mamba.readthedocs.io/) works the same
way and is faster:

```bash
mamba install -c bioconda vardictcpp
```

## Galaxy

vardictcpp is available as a Galaxy tool from the
[Galaxy Tool Shed](https://toolshed.g2.bx.psu.edu/), published under the author's own repository
(owner `mhh-hematology`). A Galaxy administrator can install it from
**Admin > Tool Management > Install and Uninstall**, by searching for `vardictcpp` in the Tool Shed.
The wrapper exposes single-sample and paired (tumor/normal) modes and the common calling options.

## Build from source

Requirements: a C++17 compiler (GCC or Clang), [CMake](https://cmake.org/) 3.15 or newer, and htslib
1.10 or newer.

```bash
git clone https://github.com/MHH-Bioinformatics-Hematology/vardictcpp.git
cd vardictcpp
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

The binary is written to `build/vardictcpp`. If htslib is not on the default search path (for example
inside a conda environment), point CMake at it:

```bash
cmake -B build -DHTSLIB_ROOT="$CONDA_PREFIX"
```

The default build targets the architecture baseline (SSE2 on x86-64, NEON on ARM), so the binary is
portable across machines of the same architecture. To tune for the build host instead, configure with
`-DVARDICTCPP_NATIVE=ON`.

## Verify the installation

```bash
vardictcpp --version
```
