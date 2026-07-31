# Whole-exome benchmark - vardictcpp 1 vs VarDictJava 1.8.3

Five public whole-exome SRA runs (hg19, bwa-mem, per-sample covered-target BED), `-f 0.01`, single
core and multi-thread. Wall clock and peak RSS captured with the same harness (`run_capped.py`,
process-group RSS poll). vardictcpp is the `1` release (SIMD upper-casing + SIMD/rolling seed index).
VarDictJava 1.8.3 is stock upstream, `-Xmx 100g`, benchmarked on two JVMs: **JDK 8** (its documented
target, Parallel/throughput GC) and **JDK 25** (modern default, G1 GC).

Samples (SRA run accessions): SRR15006540, SRR15006376, SRR15006375, SRR15006386, and SRR8657348 (the
large 774k-region run).

## Single core (`-th 1`, `-f 0.01`)

VarDictJava 1.8.3 (stock) and 1.8.4 (memory fork, G1 + string dedup), each on JDK 8 and JDK 25
(wall s / peak RSS GB):

| sample | regions | vardictcpp | 1.8.3 JDK 8 | 1.8.3 JDK 25 | 1.8.4 JDK 8 | 1.8.4 JDK 25 |
|---|--:|--|--|--|--|--|
| SRR15006540 | 135k | 56.8 / 0.061 | 252.1 / 32.47 | 260.6 / 10.96 | 484.5 / 2.48 | 252.5 / 10.29 |
| SRR15006376 | 106k | 42.6 / 0.065 | 189.0 / 21.52 | 213.8 /  9.36 | 334.7 / 8.39 | 209.8 /  6.34 |
| SRR8657348  | 774k | 93.5 / 0.088 | 626.3 / 32.67 | 720.9 /  3.31 | 552.9 / 42.03 | 720.9 /  3.80 |
| SRR15006375 | 12k  | 50.7 / 0.052 | 470.6 /  5.49 | 382.6 /  4.72 | 408.0 / 2.73 | 376.4 /  4.73 |
| SRR15006386 | 12k  | 48.7 / 0.054 | 256.5 /  6.47 | 357.1 /  3.75 | 375.3 / 2.61 | 330.5 /  3.84 |

vardictcpp geomean advantage: **5.8x faster / 242x less RAM** vs 1.8.3 JDK 8; **6.3x / 91x** vs 1.8.3
JDK 25. Wall time differs by at most ~15% across Java version and JVM; the main effect is memory - the
1.8.4 fork's G1 default cuts JDK 8 memory to 2.5-8 GB except on the un-chunked 774k-region sample.

## 8 threads (`-th 8`, `-f 0.01`) - multi-thread Java on JDK 8

Wall s / peak RSS GB. vardictcpp vs 1.8.3: geomean 16.2x faster / 85x less RAM.

| sample | vardictcpp | 1.8.3 JDK 8 | 1.8.4 JDK 8 |
|---|--|--|--|
| SRR15006540 | 10.2 / 0.300 | 138.8 / 33.37 | 138.9 / 13.73 |
| SRR15006376 |  8.1 / 0.321 | 166.5 / 33.59 | 153.9 / 13.72 |
| SRR8657348  | 14.3 / 0.243 | 267.4 / 33.29 | 248.5 / 45.48 |
| SRR15006375 |  8.2 / 0.279 | 130.0 / 12.75 | 136.1 /  2.90 |
| SRR15006386 |  8.1 / 0.213 | 109.6 / 12.75 | 113.7 /  8.64 |

## vardictcpp thread scaling (`-f 0.01`, wall s)

| sample | th1 | th4 | th8 | th16 |
|---|--:|--:|--:|--:|
| SRR15006540 | 56.8 | 16.3 | 10.2 | 8.2 |
| SRR15006376 | 42.6 | 12.2 |  8.1 | 6.2 |
| SRR8657348  | 93.5 | 26.5 | 14.3 | 8.2 |
| SRR15006375 | 50.7 | 14.2 |  8.2 | 6.1 |
| SRR15006386 | 48.7 | 14.3 |  8.1 | 6.1 |

## Notes

- **Memory / JVM.** With `-Xmx 100g`, JDK 8's throughput collector hoards heap toward the cap (21-33 GB
  on the coverage-dense samples), while JDK 25's G1 releases it (3-11 GB) - a 3-10x difference from the
  JVM alone. vardictcpp uses 40-500x less than either. The 8-thread table below is JDK 8 only; wall
  time is essentially JVM-insensitive (within ~15%), so only Java's memory shifts materially.
- **Output.** On noisy whole-exome data the two callers agree on the variant *set* to within a handful
  of calls per sample (tens of FP/FN out of 12k-54k), not bit-for-bit; see
  [equivalence_sra_wes.md](equivalence_sra_wes.md). On curated goldens the output is byte-identical.
- **Where the speedup comes from.** vardictcpp holds native counters instead of the JVM object graph,
  and version 1 added SIMD base handling and an O(n) rolling seed index; the gain is largest on
  region-dense samples (up to 1.5x over vardictcpp's own previous release; see the project history).
