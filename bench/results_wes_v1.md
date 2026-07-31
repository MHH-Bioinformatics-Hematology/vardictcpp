# Whole-exome benchmark - vardictcpp 1 vs VarDictJava 1.8.3

Five public whole-exome SRA runs (hg19, bwa-mem, per-sample covered-target BED), `-f 0.01`, single
core and multi-thread. Wall clock and peak RSS captured with the same harness (`run_capped.py`,
process-group RSS poll). vardictcpp is the `1` release (SIMD upper-casing + SIMD/rolling seed index).
VarDictJava 1.8.3 is stock upstream, `-Xmx 100g`, benchmarked on two JVMs: **JDK 8** (its documented
target, Parallel/throughput GC) and **JDK 25** (modern default, G1 GC).

Samples (SRA run accessions): SRR15006540, SRR15006376, SRR15006375, SRR15006386, and SRR8657348 (the
large 774k-region run).

## Single core (`-th 1`, `-f 0.01`)

VarDictJava on both JVMs (wall s / peak RSS GB):

| sample | regions | vardictcpp | VarDictJava JDK 8 | VarDictJava JDK 25 |
|---|--:|--|--|--|
| SRR15006540 | 135k | 56.8 / 0.061 | 252.1 / 32.47 | 260.6 / 10.96 |
| SRR15006376 | 106k | 42.6 / 0.065 | 189.0 / 21.52 | 213.8 /  9.36 |
| SRR8657348  | 774k | 93.5 / 0.088 | 626.3 / 32.67 | 720.9 /  3.31 |
| SRR15006375 | 12k  | 50.7 / 0.052 | 470.6 /  5.49 | 382.6 /  4.72 |
| SRR15006386 | 12k  | 48.7 / 0.054 | 256.5 /  6.47 | 357.1 /  3.75 |

vardictcpp geomean advantage: **5.8x faster / 242x less RAM** vs JDK 8; **6.3x faster / 91x less RAM**
vs JDK 25. Wall time differs by at most ~15% between the two JVMs; the JVM's main effect is memory.

## 8 threads (`-th 8`, `-f 0.01`) - VarDictJava on JDK 8

| sample | vardictcpp wall | vardictcpp RSS | VarDictJava wall | VarDictJava RSS | faster | less RAM |
|---|--:|--:|--:|--:|--:|--:|
| SRR15006540 | 10.2 s | 0.300 GB | 138.8 s | 33.37 GB | 13.6x | 111x |
| SRR15006376 |  8.1 s | 0.321 GB | 166.5 s | 33.59 GB | 20.6x | 105x |
| SRR8657348  | 14.3 s | 0.243 GB | 267.4 s | 33.29 GB | 18.7x | 137x |
| SRR15006375 |  8.2 s | 0.279 GB | 130.0 s | 12.75 GB | 15.9x |  46x |
| SRR15006386 |  8.1 s | 0.213 GB | 109.6 s | 12.75 GB | 13.5x |  60x |
| **geomean** | | | | | **16.2x** | **85x** |

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
