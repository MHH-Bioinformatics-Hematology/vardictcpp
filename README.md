# vardictcpp

A C++17 port of [VarDict](https://github.com/AstraZeneca-NGS/VarDictJava) (AstraZeneca-NGS),
built with htslib. Goal: a memory-lean, fast native implementation of the VarDict amplicon/somatic
variant caller. This repository is the **staged port**; see *Parity status* for what is implemented
today.

## Why

VarDictJava can use very large amounts of RAM (reported ~200 GB with `-f 0`). Two causes, confirmed
by profiling the Java version:

1. **JVM heap over-allocation** — the live set is far smaller than resident memory; with a large
   `-Xmx` the throughput collector lets garbage accumulate, so RSS tracks `-Xmx`.
2. **Per-region live-set scaling** — VarDict builds per-position variation maps for the whole region
   at once, so a whole-chromosome `-R` interval holds a genome-scale object graph.

A native C++ implementation removes (1) entirely (no managed heap; memory is released as soon as a
region's structures go out of scope) and lets us bound (2) by streaming and region windowing.

## Build

Requires a C++17 compiler, CMake ≥ 3.15, and htslib (headers + lib). With htslib from conda:

```bash
cmake -S . -B build -DHTSLIB_ROOT=$CONDA_PREFIX
cmake --build build -j
```

`build/vardictcpp` links htslib via rpath, so it runs without `LD_LIBRARY_PATH`.

## Usage

```bash
vardictcpp -G ref.fa -b in.bam -N sample -R chr:start-end -f 0.01
vardictcpp -G ref.fa -b in.bam -N sample -c 1 -S 2 -E 3 -g 4 panel.bed -f 0
vardictcpp -G ref.fa -b in.bam -N sample -R chr:1-1000000 --chunk 50000   # bound memory
```

Output is VarDict simple-mode TSV (same 36-column order as `SimpleOutputVariant`). `--chunk N`
splits regions longer than `N` bp into consecutive windows, bounding peak memory independently of
interval length (identical mechanism to the `--chunk` flag added to VarDictJava in this project).

## Architecture (mirrors the Java package layout)

| C++ file | Ports from (VarDictJava) | Purpose |
|---|---|---|
| `src/config.hpp` | `Configuration.java` | run options |
| `src/region.hpp` | `RegionBuilder.java` | region model + `splitLongRegions` (`--chunk`) |
| `src/reference.{hpp,cpp}` | `data/ReferenceResource.java` | indexed FASTA via htslib `faidx` |
| `src/variation.hpp` | `variations/Variation.java` | per-allele accumulator |
| `src/cigar_parser.{hpp,cpp}` | `modules/CigarParser.java` | BAM read → CIGAR walk → per-position counts (**streaming**, one record at a time) |
| `src/tovars.{hpp,cpp}` | `modules/ToVarsBuilder.java` + `postprocessmodules/SimplePostProcessModule.java` | per-position variant calling + metrics |
| `src/printer.{hpp,cpp}` | `printers/SimpleOutputVariant.java` | TSV output |

## Parity status

**Implemented (validated):** single-sample counting core — BAM iteration with SAM-flag / mapq
filtering and optional duplicate removal; CIGAR M/=/X/I/D/N/S handling; per-position ref coverage and
per-allele counts (strand, base-quality, mapping-quality, NM, hi/lo-quality, pstd/qstd) using
VarDict's exact read-position convention (distance to nearest read end); `-f` frequency filter;
**`isGoodVar`** quality gate for default (non-pileup) mode; strand-bias flag; AF/PMean/PStd/QMean/QStd/
MQ/Sig_Noise/HiAF; **`findMSI`** (MSI / MSI_NT / shift3) for the SNV/MNP path; 20 bp reference flanks;
the exact 36-column output formatting (`value==0 ? "0"` rules); `--chunk` windowing; streaming
per-region memory release.

**Validation vs VarDictJava 1.8.3:**
- *Pileup counting* (1 Mb / 300× synthetic, matched SNV alleles): Depth exact **97.7 %** (mean abs
  diff 0.03 reads), AltDepth exact **93.3 %** (0.14 reads).
- *Default simple mode* (hg19 panel, `-f 0.01`): **6 of VarDictJava's 9 calls reproduce byte-for-byte
  across all 36 columns** (counts, AF, PMean/QMean, MSI, flanks, Seg, Duprate). The 3 that differ are
  an MNP (`TCC>ACG`, needs `adjustMNP`), an insertion needing realignment left-normalization, and one
  SNV off by a single alt read (see below).

**Not yet ported (needed for full byte-for-byte parity):**

- The ±1-read residual on some positions: VarDict's base-quality / neighbour handling in the matching
  increment (MNV growth, insertion-edge count adjustment). (`skipOverlappingReads` is *off* by default
  — only under `-u`/`-UN` — so it is not the cause.)
- **`adjustMNP`** — merges adjacent SNVs into MNPs (the `TCC>ACG` case).
- Local **realignment** (`VariationRealigner`: realigndel/ins/lgdel/lgins) and soft-clip consensus —
  required for correct indel calls and left-normalization.
- **Structural variants** (`StructuralVariantsProcessor`).
- **Somatic** (paired) and **amplicon** modes; `--fisher` (hypergeometric + Brent `zeroin`).

These map 1:1 onto the remaining Java modules (see the table) and are the roadmap to full parity.

## License

Mirrors the upstream VarDict license (MIT). This is an independent reimplementation for research use.
