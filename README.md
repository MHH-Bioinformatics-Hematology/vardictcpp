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
vardictcpp -G ref.fa -b in.bam -N sample -c 1 -S 2 -E 3 -g 4 panel.bed -f 0 --th 8
vardictcpp -G ref.fa -b in.bam -N sample -R chr:1-1000000 --chunk 50000   # bound memory
```

Output is VarDict simple-mode TSV (same 36-column order as `SimpleOutputVariant`). `--chunk N`
splits regions longer than `N` bp into consecutive windows, bounding peak memory independently of
interval length (identical mechanism to the `--chunk` flag added to VarDictJava in this project).
`--th N` (alias `--threads`) processes regions across `N` worker threads with ordered streaming
output (mirrors VarDictJava's parallel mode); output is bit-identical to single-threaded. On the
698-region hg19 panel at 8 threads: **0.52 s / 0.05 GB** vs VarDictJava **1.88 s / 1.47 GB**.

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

Includes VarDict's **MNV growth** (adjacent mismatches, bridging up to `vext` matching bases, grown
into one `&`-joined description string, e.g. `A&CG` → `TCC>ACG` Complex) and the MNP `mnp` map.

Now also ports: **soft-clip processing** (mis-softclip re-matching that converts reference-matching
clipped bases to coverage, plus per-position consensus in `softClips5End`/`softClips3End` for
realignment) and **`adjustMNP`** (merges partial SNVs into the MNP they belong to via `adjCnt` and
removes them). ExtraAF is derived from `extracnt`.

**Validation vs VarDictJava 1.8.3:**
- *Pileup counting* (1 Mb / 300× synthetic, matched SNV alleles): Depth exact **97.9 %** (mean abs
  diff 0.027 reads), AltDepth exact **99.7 %** (0.006 reads) — MNV growth brought AltDepth up from
  93.3 %.
- *Default simple mode* (hg19 panel, `-f 0.01`): **7 of VarDictJava's 9 calls reproduce byte-for-byte
  across all 36 columns**; the single-sample **SNV path is effectively at parity** (the last ±1-read
  residual was resolved by soft-clip re-matching). The `TCC>ACG` MNP now has correct Ref/Alt/type/MSI
  and AltDepth (Depth still differs — see coverage note). The port currently emits extra indel/edge
  calls that VarDict's realignment reassigns or removes (below).

Also ports the **small-indel realignment engine** (`realigner.cpp`): `realignins`/`realigndel` with
`findMM3`/`findMM5`, `findconseq`, `ismatch`, `joinRef`, `adjCnt`/`adjRefCnt`/`adjRefFactor`, and
insertion left-normalization (`adjInsPos`). These attribute nearby mismatch SNVs and soft-clip
consensus reads to an indel and merge duplicate representations. They run in VarDict's order
(`realigndel` → `realignins` → `adjustMNP`) and are functioning + non-regressing.

**The complete soft-clip realignment engine is ported and enabled**, in VarDict's order:
`adjustMNP` → `realigndel` → `realignins` → `realignlgdel` → `realignlgins30` → `realignlgins`.
This includes `findbp`, `findbi`, `find35match`, the reference **k-mer seed index** + `findMatch`,
`findconseq` (with the `B_A7`/`B_T7` poly-A/T guard), and the soft-clip-consensus / mismatch
reassignment that removes indel-explained SNVs. Reference is loaded with VarDict's ±1200 window.

**Not yet ported — `CigarModifier` (the precise cause of the remaining false-positives):**

Diagnosed to ground truth: the remaining panel FPs are reads whose CIGAR VarDict **rewrites before
counting** and this port does not. E.g. at `chr3:47538004` the reads are `64S36M112S` (a 36 bp mapped
island, mate unmapped) in a **poly-T homopolymer**; the 36 M has a single G>C. `modules/CigarModifier`
(787 lines: leading/trailing D/I normalization, chimeric-seed clip removal, `captureMisSoftlyMS`/
`captureMisSoftly3Mismatches`, `combineDigSDigM`/`combineBeginDigM`, and the indel-collapse loop)
reshapes such reads so they never produce the SNV. VarDict pileup counts **nothing** there; this port
counts `G>C 29/29`. `CigarModifier` runs on *every* read at the top of `parseCigar`, so it must be
ported carefully (a defect changes all counting, not just these rows) — it is deliberately left for a
dedicated pass rather than risking the verified pipeline.

Also remaining: the **discordant/chimeric SV subsystem** (`StructuralVariantsProcessor` +
`SVStructures` clustering) for SV *output*; faithful **distributed coverage** at indel/MNP-dense
positions (the `TCC>ACG` Depth 46 vs 35); **somatic** (paired) and **amplicon** modes; `--fisher`.

Ported & enabled: CIGAR parse → MNV/MNP → soft-clip → **full small + large indel realignment** →
call/format. Remaining: the discordant/chimeric SV subsystem, distributed coverage, somatic/amplicon,
fisher.

## License

Mirrors the upstream VarDict license (MIT). This is an independent reimplementation for research use.
