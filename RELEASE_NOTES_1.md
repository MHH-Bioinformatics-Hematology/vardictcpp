# vardictcpp 1

First release of **vardictcpp**, a C++17 reimplementation of
[VarDict](https://github.com/AstraZeneca-NGS/VarDictJava) (AstraZeneca-NGS VarDictJava 1.8.3) built on
htslib. It reproduces VarDict's variant calls byte-for-byte while running several times faster and using
roughly an order of magnitude less memory.

## Highlights
- **Byte-identical to VarDictJava 1.8.3** across all ported modes on five whole-exome samples
  (91,288 variant rows, 0 FP / 0 FN) plus per-mode golden fixtures enforced in CI.
- **Full mode coverage**: simple, amplicon, `--fisher`, structural variants (`<INV>`, `<DEL>`, `<DUP>`
  including discordant pairs and inter-chromosomal fusion clusters), paired somatic (both BAMs compared,
  including the merged-region `combineAnalysis` refinement), and splice-junction handling.
- **Independently accuracy-validated** against a public truth set (GIAB HG002, NIST v4.2.1), where
  vardictcpp and VarDictJava produce identical calls.
- **Faster and far leaner** than stock VarDictJava (numbers below).
- **Streaming and region-parallel** (`-th`) with output identical to single-threaded; optional
  `--chunk` bounds memory on very large regions.
- **Portable SIMD** (SSE2 on x86-64, NEON on ARM/Apple Silicon, scalar fallback), every kernel
  bit-for-bit equal to its scalar form. CI builds and tests on Linux and macOS, x86-64 and ARM.
- Drop-in CLI: accepts VarDict's option syntax and emits VarDict's TSV column order (pipe straight into
  `teststrandbias.R` / `var2vcf_valid.pl`).

## Correctness (parity with VarDictJava 1.8.3)
- **Byte-identical on all five WES samples tested**, zero differing lines versus single-threaded
  VarDictJava 1.8.3 (multi-threaded Java output is non-deterministic, so single-threaded is the
  reference):

  | Sample | Target regions | Variant rows | Byte-identical |
  |---|---:|---:|---:|
  | SRR15006386 | 11,732 | 4,014 | 4,014 (100%) |
  | SRR15006375 | 12,292 | 4,056 | 4,056 (100%) |
  | SRR15006376 | 105,949 | 12,739 | 12,739 (100%) |
  | SRR15006540 | 135,145 | 16,031 | 16,031 (100%) |
  | SRR8657348 (MV4-11 / CCLE) | 773,667 | 54,448 | 54,448 (100%) |
  | **Total** | | **91,288** | **91,288 (100%)** |

- Reaching this drove the full structural-variant subsystems, large-indel coverage reloads, the
  reference `SEED_1` extent truncation, and the CigarModifier soft-clip/homopolymer fixes.
- Every divergence found during development was triaged read-by-read against **both** VarDictJava and
  the original Perl: each proved to be a vardictcpp bug corrected toward the reference, so
  `docs/DIVERGENCES.md` has no open entries. Two latent correctness bugs surfaced in the process (an
  output-buffer overflow that silently dropped large-allele variant lines, and a soft-clip read filter
  applied to the pre-modification CIGAR) and are fixed.
- Other WES samples are not exhaustively verified; new data may surface further edge cases.

## Accuracy against ground truth
Validated on the GIAB HG002 chromosome-20 exome (Agilent SureSelect v5) against the NIST v4.2.1
benchmark with `rtg vcfeval` over callable confident regions (germline operating point: PASS,
allele frequency >= 0.2). vardictcpp and VarDictJava emit **byte-identical VCFs** here as well
(8,244 variants, 0 differing lines), so both share exactly the same accuracy:

| Variant class | Precision | Recall | F1 |
|---|---:|---:|---:|
| SNV | 0.996 | 0.974 | 0.985 |
| Indel | 0.877 | 0.734 | 0.799 |
| All | 0.984 | 0.966 | 0.975 |

## Performance vs stock VarDictJava 1.8.3
Real whole-exome data (public SRA runs aligned to hg19, covered-target BED, `-f 0.01`), CPU-pinned,
median of 5 clean replicate runs, peak RSS via `/usr/bin/time -v`, VarDictJava on JDK 25 with `-Xmx 8g`:

| Threads | Speedup (Java / cpp) | Peak RSS (cpp) | Peak RSS (Java) | RSS reduction |
|---|---|---|---|---|
| 1 | 3.7 to 4.8x | 74 to 122 MB | ~1.0 to 1.3 GB | 8 to 17x |
| 8 | 9 to 11x | 252 to 358 MB | ~2.1 to 2.6 GB | 7 to 8x |

The memory advantage grows with region size: at `-f 0` VarDict retains every covered base times allele,
whereas vardictcpp streams one record at a time and stores the reference as disjoint windows matching
Java's reference map, so a distant realignment breakpoint never gap-fills a multi-megabase span.

## Build & test
- Requires a C++17 compiler, CMake >= 3.15, and htslib:
  `cmake -S . -B build -DHTSLIB_ROOT=$CONDA_PREFIX && cmake --build build -j`
- `bash test/run_tests.sh`: parity against self-contained VarDictJava golden fixtures (all modes).
- GitHub Actions CI (Linux + macOS, x86-64 + ARM): gcc/clang build + parity tests, plus a runtime and
  memory benchmark vs VarDictJava.
- `-DVARDICTCPP_NATIVE=ON` enables `-march=native`; the default build stays portable.

## License
MIT, retaining the upstream AstraZeneca-NGS VarDictJava copyright. See `LICENSE` and `NOTICE`.

## Known limitations
- Real-WES parity is byte-identical on the five verified samples; other samples are not exhaustively
  verified.
- The ground-truth accuracy benchmark is chromosome 20 only (to bound compute); it characterises
  VarDict's accuracy, which the port inherits byte-identically, and is not a novel-accuracy claim.
- `--chunk` is opt-in; without it a single very large `-R` region is one work unit (as in VarDictJava).
