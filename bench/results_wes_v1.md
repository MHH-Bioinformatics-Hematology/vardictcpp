# Whole-exome benchmark - vardictcpp 1 vs VarDictJava 1.8.3

Five public whole-exome SRA runs (hg19, bwa-mem, per-sample covered-target BED), `-f 0.01`, single
core and multi-thread. Wall clock and peak RSS captured with the same harness (`run_capped.py`,
process-group RSS poll). vardictcpp is the `1` release (SIMD upper-casing + SIMD/rolling seed index).
VarDictJava 1.8.3 is stock upstream on JDK 8 (its documented target JVM), `-Xmx 100g`.

Samples (SRA run accessions): SRR15006540, SRR15006376, SRR15006375, SRR15006386, and SRR8657348 (the
large 774k-region run).

## Single core (`-th 1`, `-f 0.01`)

| sample | regions | vardictcpp wall | vardictcpp RSS | VarDictJava wall | VarDictJava RSS | faster | less RAM |
|---|--:|--:|--:|--:|--:|--:|--:|
| SRR15006540 | 135k | 56.8 s | 0.061 GB | 252.1 s | 32.47 GB | 4.4x | 532x |
| SRR15006376 | 106k | 42.6 s | 0.065 GB | 189.0 s | 21.52 GB | 4.4x | 331x |
| SRR8657348  | 774k | 93.5 s | 0.088 GB | 626.3 s | 32.67 GB | 6.7x | 371x |
| SRR15006375 | 12k  | 50.7 s | 0.052 GB | 470.6 s |  5.49 GB | 9.3x | 106x |
| SRR15006386 | 12k  | 48.7 s | 0.054 GB | 256.5 s |  6.47 GB | 5.3x | 120x |
| **geomean** | | | | | | **5.8x** | **242x** |

## 8 threads (`-th 8`, `-f 0.01`)

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

- **Memory / JVM.** The VarDictJava peak RSS above is JDK 8 (VarDict's target JVM) with `-Xmx 100g`,
  so it includes GC head-room. A modern JVM (G1, e.g. JDK 21+) lowers VarDictJava's single-core peak
  RSS to roughly 3-11 GB - still 40-180x above vardictcpp. Wall time is essentially JVM-insensitive.
- **Output.** On noisy whole-exome data the two callers agree on the variant *set* to within a handful
  of calls per sample (tens of FP/FN out of 12k-54k), not bit-for-bit; see
  [equivalence_sra_wes.md](equivalence_sra_wes.md). On curated goldens the output is byte-identical.
- **Where the speedup comes from.** vardictcpp holds native counters instead of the JVM object graph,
  and version 1 added SIMD base handling and an O(n) rolling seed index; the gain is largest on
  region-dense samples (up to 1.5x over vardictcpp's own previous release; see the project history).
