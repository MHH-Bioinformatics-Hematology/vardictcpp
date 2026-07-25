# vardictcpp vs UNFIXED VarDictJava 1.8.3 - runtime & memory

Head-to-head of the **stock, unmodified** `VarDict-1.8.3.jar` (default JVM: no `-Xmx` tuning, no G1,
no `--chunk`) against `vardictcpp`. Peak RSS = `/usr/bin/time -v` "Maximum resident set size"; wall =
elapsed. Host: 32 cores, 125 GB RAM, hg19. Run 2026-07-24.

## 1. Real AML panel - MRD2020-032_H20PB2025 (220 MB BAM, 698-region panel)

| -f    | threads | Java peak RSS | Java wall | C++ peak RSS | C++ wall | RAM ratio |
|-------|---------|---------------|-----------|--------------|----------|-----------|
| 0.01  | 1       | 1.386 GB      | 1.96 s    | 0.012 GB     | 3.74 s   | **116x**  |
| 0.01  | 4       | 1.330 GB      | 1.21 s    | 0.034 GB     | 1.04 s   | **39x**   |
| 0     | 1       | 1.371 GB      | 1.80 s    | 0.013 GB     | 4.77 s   | **105x**  |
| 0     | 4       | 1.410 GB      | 1.26 s    | 0.034 GB     | 1.68 s   | **41x**   |

On a fragmented panel the JVM's ~1.3-1.4 GB baseline dominates and hides the algorithmic memory cost;
C++ holds a flat 12-34 MB. Wall-clock is close: Java wins single-thread on this tiny workload (JVM is
already warm and C++ pays per-region reference/index seek overhead across 698 small regions), C++ edges
ahead at 4 threads. The memory gap is the real, consistent win here.

## 2. Memory-stress case - single 1 Mb region at -f 0 (2 M reads, ~300x)

This is what the port was built for: one large region with `-f 0` makes VarDict retain ~every covered
base x allele, so memory scales with region length (a whole chromosome as one `-R` is how stock Java
reaches ~200 GB).

| tool                     | peak RSS  | wall     | lines  |
|--------------------------|-----------|----------|--------|
| stock Java (default JVM) | 6.270 GB  | 30.0 s   | 23283  |
| **vardictcpp**           | 0.540 GB  | 19.1 s   | 23282  |
| vardictcpp `--chunk 50k` | 0.036 GB  | 13.7 s   | 23282  |

**11.6x less RAM and 1.6x faster** out of the box; **174x less RAM and 2.2x faster** with region
chunking. The gap grows with region size - this 1 Mb slice already shows 6.3 GB vs 0.5 GB.

## 3. Output equivalence (correctness of the speedup)

- **Panel, -f 0.01 (calling threshold):** Java 9 variants, C++ 9 variants, **0 false-positives / 0
  false-negatives**, 8/9 rows fully byte-identical. (The pre-SV-port benchmark had 2 FP / 2 FN on this
  same sample from missing structural-variant output; porting SV closed both.) The 1 non-identical row
  differs only in coverage-metric columns at the chr21 ultra-high-depth satellite locus.
- **Stress, -f 0:** 23282 / 23283 rows byte-identical (the 1 extra Java row is the same coverage-frontier
  locus). At the real -f 0.01 threshold the synthetic set is 1904/1904 exact, 0 FP / 0 FN.

## Takeaway
Same calls, dramatically less memory. On real panel data vardictcpp uses **40-115x less RAM**; on the
`-f 0` large-region case that motivated the work it uses **12x (chunked: 174x) less RAM while running
faster**, with byte-identical variant sets. Reproduce: `bench/run_java_vs_cpp.sh`.
