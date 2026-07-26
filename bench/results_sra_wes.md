# UNFIXED VarDictJava 1.8.3 vs vardictcpp on 3 SRA DNA runs (hg19, bwa mem, covered-target BED, -f 0.01)
Generated 2026-07-25 10:35:56. Peak RSS via /usr/bin/time -v "Maximum resident set size".

| sample | type | threads | tool | peak RSS | wall | vars |
|---|---|---|---|---|---|---|
| SRR15006540 | DNA | 1 | stock Java | 1.886 GB | 4:01.94 | 16031 |
| SRR15006540 | DNA | 1 | vardictcpp | 0.062 GB | 0:58.32 | 16026 |
| SRR15006540 | DNA | 4 | stock Java | 2.768 GB | 1:58.07 | 16031 |
| SRR15006540 | DNA | 4 | vardictcpp | 0.191 GB | 0:15.61 | 16026 |
| SRR15006376 | DNA | 1 | stock Java | 5.012 GB | 3:01.88 | 12739 |
| SRR15006376 | DNA | 1 | vardictcpp | 0.073 GB | 0:43.27 | 12733 |
| SRR15006376 | DNA | 4 | stock Java | 4.153 GB | 1:39.99 | 12739 |
| SRR15006376 | DNA | 4 | vardictcpp | 0.198 GB | 0:11.95 | 12733 |
| SRR8657348 | DNA | 1 | stock Java | 4.171 GB | 8:22.64 | 54448 |
| SRR8657348 | DNA | 1 | vardictcpp | 0.119 GB | 1:43.17 | 54454 |
| SRR8657348 | DNA | 4 | stock Java | 5.396 GB | 3:11.28 | 54448 |
| SRR8657348 | DNA | 4 | vardictcpp | 0.175 GB | 0:27.67 | 54454 |

## Summary (final optimized vardictcpp b34f448 vs unfixed VarDictJava 1.8.3)

**Speed** (wall-clock, speedup = Java/C++): mean **4.4x faster single-thread**, **7.6x faster at 4 threads**.
Per sample th4: SRR15006540 7.6x, SRR15006376 8.4x, SRR8657348 6.9x. C++ is faster than Java at BOTH
thread counts on all 3 samples (before the optimization work C++ was 2-14x SLOWER on this same WES data).

**Memory** (peak RSS, reduction = Java/C++): **14-69x less** (C++ 0.06-0.20 GB vs Java 1.9-5.4 GB).

**Correctness** (-f 0.01, -th 4): ~95% of variant rows byte-identical; small residual per sample
(SRR15006540 5 FP/10 FN, SRR15006376 1 FP/7 FN, SRR8657348 43 FP/37 FN of ~16k-54k calls) in the known
edge-case classes (SV/large-indel representation, distributed indel coverage, homopolymer/MNV). The loop
drove these down from the pre-optimization 214/163/595 FP earlier; they are not yet zero on full WES.

Figures: bench/figures/fig_speed.{pdf,svg,png}, fig_memory.{pdf,svg,png} (make_figures.py).
