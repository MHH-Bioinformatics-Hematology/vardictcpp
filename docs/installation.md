# Installation

vardictcpp runs on Linux (x86-64 and ARM) and macOS on Apple Silicon (ARM). Its only runtime
dependency is [htslib](https://github.com/samtools/htslib).

!!! note "Intel (x86-64) Macs"
    vardictcpp is architecture-portable, and we believe it also runs on Intel-based Macs; we try to
    provide a conda package for that platform. It is not officially supported, however: Apple has
    announced that macOS Tahoe (26) is the last release to support Intel Macs, with macOS 27 requiring
    Apple Silicon, so we do not support Intel macOS officially. Do not be surprised if a newer
    vardictcpp release stops working on an Intel Mac, and note that Intel-Mac-specific breakage will
    not be fixed.

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

vardictcpp is available as a Galaxy tool. The wrapper is maintained at
[MHH-Bioinformatics-Hematology/galaxytools](https://github.com/MHH-Bioinformatics-Hematology/galaxytools)
and distributed through the author's [Galaxy Tool Shed](https://toolshed.g2.bx.psu.edu/) channel
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
