# vardictcpp

[![CI](https://github.com/MHH-Bioinformatics-Hematology/vardictcpp/actions/workflows/ci.yml/badge.svg)](https://github.com/MHH-Bioinformatics-Hematology/vardictcpp/actions/workflows/ci.yml)

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

### Portability and SIMD

The default build is **portable**: no `-march=native`, so the binary runs on any CPU of its
architecture and builds on Linux and macOS, x86-64 and ARM / Apple Silicon (arm64). Hot loops use
explicit SIMD through a small portable layer (`src/simd.hpp`) that selects a backend from the
architecture *baseline* instruction set, so no ISA flag is needed and the result is bit-identical to
the scalar path:

- **x86-64** → SSE2 (guaranteed on every x86-64 CPU)
- **ARM / Apple Silicon** → NEON
- anything else → a scalar fallback

`vardictcpp --version` prints the version and the compiled-in backend (e.g.
`vardictcpp 1 (SIMD backend: NEON)`). For a locally built, non-distributed binary you can add
host-specific tuning with `-DVARDICTCPP_NATIVE=ON` (adds `-march=native`, or `-mcpu=native` on ARM);
leave it off for anything you ship, and never use it for a Bioconda build.

CI builds and runs the parity suite on the full matrix (Linux and macOS, x86-64 and arm64, gcc and
clang), so portability across those targets is enforced on every push.

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
output (mirrors VarDictJava's parallel mode); output is bit-identical to single-threaded.

## Whole-exome benchmark

Five public whole-exome SRA runs (hg19, bwa-mem, per-sample covered-target BED), `-f 0.01`, vs stock
VarDictJava 1.8.3. Wall clock and peak RSS from the same harness; full numbers and method in
[bench/results_wes_v1.md](bench/results_wes_v1.md).

**Single core (`-th 1`)**

| sample | regions | vardictcpp | VarDictJava 1.8.3 | faster | less RAM |
|---|--:|--|--|--:|--:|
| SRR15006540 | 135k | 56.8 s / 0.06 GB | 252.1 s / 32.5 GB | 4.4x | 532x |
| SRR15006376 | 106k | 42.6 s / 0.07 GB | 189.0 s / 21.5 GB | 4.4x | 331x |
| SRR8657348  | 774k | 93.5 s / 0.09 GB | 626.3 s / 32.7 GB | 6.7x | 371x |
| SRR15006375 | 12k  | 50.7 s / 0.05 GB | 470.6 s /  5.5 GB | 9.3x | 106x |
| SRR15006386 | 12k  | 48.7 s / 0.05 GB | 256.5 s /  6.5 GB | 5.3x | 120x |
| **geomean** | | | | **5.8x** | **242x** |

**8 threads (`-th 8`)**

| sample | vardictcpp | VarDictJava 1.8.3 | faster | less RAM |
|---|--|--|--:|--:|
| SRR15006540 | 10.2 s / 0.30 GB | 138.8 s / 33.4 GB | 13.6x | 111x |
| SRR15006376 |  8.1 s / 0.32 GB | 166.5 s / 33.6 GB | 20.6x | 105x |
| SRR8657348  | 14.3 s / 0.24 GB | 267.4 s / 33.3 GB | 18.7x | 137x |
| SRR15006375 |  8.2 s / 0.28 GB | 130.0 s / 12.8 GB | 15.9x |  46x |
| SRR15006386 |  8.1 s / 0.21 GB | 109.6 s / 12.7 GB | 13.5x |  60x |
| **geomean** | | | **16.2x** | **85x** |

VarDictJava here is JDK 8 (its target JVM) with `-Xmx 100g`, so its peak RSS includes GC head-room; a
modern JVM (G1) lowers it to ~3-11 GB single-core, still 40-180x above vardictcpp, and wall time is
JVM-insensitive. On this noisy WES data the two callers agree on the variant set to within a handful of
calls per sample (see [bench/equivalence_sra_wes.md](bench/equivalence_sra_wes.md)); on curated goldens
output is byte-identical. For reference, on a small 698-region panel at 8 threads vardictcpp runs in
**0.52 s / 0.05 GB** vs VarDictJava **1.88 s / 1.47 GB**.

## CLI compatibility

vardictcpp accepts **VarDictJava 1.8.3's complete option set** (62 options) with the same
commons-cli syntax, including single-dash multi-char options (`-th 8`, `-VS STRICT`, `-DP`, `-mfreq`).
Options that drive the single-sample pipeline are acted on; the rest are parsed (VarDict-compatible)
even where not yet wired. Unsupported *modes* refuse cleanly instead of mis-calling: `-a/--amplicon`,
`--fisher`, and somatic (two BAMs `-b 't|n'`).

## Testing

- **vardictcpp:** `bash test/run_tests.sh` — runs against a self-contained fixture (20 kb reference +
  ~26k simulated reads, `test/data/`) and asserts parity with a golden output from stock VarDictJava:
  variant set identical (0 FP/FN), Depth/AltDepth exact, byte-identity ≥ 95%. Currently **PASS** (100%
  byte-identical on the fixture).
- **VarDictJava (the memory branch):** all **92** integration test cases pass — run via the standalone
  `IntegrationRunner` (the gradle/TestNG path is offline-blocked here), confirming the memory changes
  are byte-identical to stock.

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

**Validation vs VarDictJava 1.8.3** (1 Mb / 300× synthetic, default simple mode, **1807 variants**):
- **Variant set: exact — 0 false-positives, 0 false-negatives.**
- **Full-row byte-identical: 98.3 %** (1776 / 1807 rows match VarDictJava across all 36 columns).
- **Depth / AltDepth: 100 %** exact (mean abs diff 0.001 / 0.002 reads).
- The remaining ~1.7 % of rows are multi-variant / MNP-adjacent edge cases (genotype1 should be the
  *dominant* variant at the position, and coverage differs by a few reads where variants overlap).

Getting there required the coupled **CigarModifier + adjSNV** pair (read-end mismatch → soft-clip →
merged back into the adjacent SNV), verified read-by-read against instrumented VarDict, plus the exact
genotype rule (genotype1 = reference allele when its frequency ≥ `-f`, else the variant).
- *Default simple mode* (hg19 panel, `-f 0.01`): the C++ output is now the **same 9-variant set as
  VarDictJava with ZERO false-positives**, of which **7 of 9 reproduce byte-for-byte across all 36
  columns**. The last false-positives were removed by matching VarDict's read filter (drop reads
  soft-clipped at both ends, leading clip 10-99 bp + trailing clip ≥10 bp, and supplementary
  alignments) — found with an instrumented-Java per-read CIGAR diff harness, which also disproved the
  earlier SV/CigarModifier hypotheses. The 2 non-byte-identical rows differ only in specific
  coverage-accounting columns: the `TCC>ACG` MNP's `Depth` (the 11 overlapping indel reads;
  `createInsertion` attributes insertion coverage to position+1) and the `T>TC` insertion's
  `RefFwd`/`RefRev`/`HiCov` (`createInsertion` + `calcHicov` reconciliation). Ref/Alt/Depth/AltDepth/
  AF/MSI/genotype on both rows match.

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

**Deletion output now matches VarDict exactly** (`ToVarsBuilder` + `proceedVrefIsDeletion`): a
deletion is anchored one base 5′ of its stored position (`refallele = ref[p-1] + deleted bases`,
`varallele = ref[p-1]`, `startPosition--`), gets its microsatellite `msi`/`shift3`/`msint` from the
deleted-unit-vs-flank comparison, the `genotype1/-N` genotype form, and `Deletion` classification via
a full `varType()` port. Per-read `NM` is `edit_distance − (I+D length)` so indel gaps aren't counted
as mismatches. On the deduplicated MRD data this makes deletion rows byte-identical to VarDict on
30 of 32 columns (only the indel-position `Depth`/`AF` distributed-coverage accounting differs) and
eliminates every prior deletion false-negative.

**Validation on UMI-deduplicated MRD data** (the caller's real production input, 12 samples,
`bench/results_dedup.md`): **variant set exact on all 12 samples (0 FP / 0 FN)** and **11 of 12
fully byte-identical** (all 36 columns). The two dedup false-positives are now closed — by the `-m`
mismatch read filter (a poly-G artifact read) and the leading soft-clip + short-match + indel
`CigarModifier` rule (a spurious deletion behind a soft-clip) — and the shared deletion is fully
byte-identical after porting the deletion reference-coverage increment plus the
`beginDigitMNumberIorDNumberM` reshaping it depends on. The **one** remaining non-identical row is an
ultra-high-coverage (~1.4 M) insertion whose `Ref{Fwd,Rev}`/`HiCov`/`Sig_Noise`/`HiAF`/`PStd`/`QStd`
come from the un-ported `createInsertion` + `calcHicov` reconciliation; its AF still matches to four
decimals. Performance vs VarDictJava: **peak RSS 39× lower on average / 91× at peak** (2.57 GB →
0.028 GB on a repeat-dense sample) and **wall-clock 13× faster on average / 16× on the slowest
sample** (230 s → 14.4 s).

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

## Input validation and error messages

vardictcpp fails fast with a clear, single-line `vardictcpp: <what>` message (and a non-zero exit)
instead of crashing or silently producing no output. It reports: a malformed `-R` region (a bare
chromosome or non-numeric bounds, with the expected `chr:start-end` form), BED problems with the line
number and reason (too few columns for `-c/-S/-E/-g`, non-integer coordinates, `start > end`), a
non-numeric value for any numeric option (naming the option), missing `-G`/`-b`/`-N`, and a missing or
unindexed reference/BAM. As a preflight it also checks every requested chromosome against the
reference and BAM, warning (or erroring, if none match) on naming mismatches such as `chr7` vs `7`.

## Acknowledgments

vardictcpp is a port of **VarDict** and **VarDictJava** by the AstraZeneca-NGS team. The algorithms,
option set and output format are theirs; please cite the original VarDict paper when you use this
software:

> Lai Z, Markovets A, Ahdesmaki M, et al. *VarDict: a novel and versatile variant caller for
> next-generation sequencing in cancer research.* Nucleic Acids Research (2016) 44(11):e108.
> doi:10.1093/nar/gkw227

Upstream projects:
- VarDictJava: https://github.com/AstraZeneca-NGS/VarDictJava
- VarDict (Perl): https://github.com/AstraZeneca-NGS/VarDict

## License

vardictcpp is released under the **MIT License**, the same license as VarDictJava. Because this is a
port whose design and output are derived from VarDictJava, the original AstraZeneca-NGS copyright and
permission notice is retained alongside the port's, as the MIT License requires. See the
[LICENSE](LICENSE) file for the full text.

- Copyright (c) 2019 AstraZeneca - NGS Team (original VarDict / VarDictJava)
- Copyright (c) 2026 Joachim Wolff and the vardictcpp contributors

Third-party components: this project links **htslib** (MIT/Expat) for BAM/CRAM and FASTA access; no
VarDictJava source or its bundled libraries (JRegex, Commons CLI, Commons Math, htsjdk) are included.
