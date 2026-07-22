# Benchmark: stock-Java vs optimized-Java vs vardictcpp

Synthetic 1 Mb region at 300x coverage (wgsim on hg19 chr22:20,000,000-21,000,000, contig renamed
`sim`), single region `-R sim:1-1000000`, single-threaded. Peak RSS = /usr/bin/time -v "Maximum
resident set size". All three tools produce the same variant *content* (Java rows are byte-identical
across stock/opt; vardictcpp emits the same calls plus a few more low-level candidates -> see
count-agreement note).

| tool / setting              | -f 0.01 RSS | -f 0.01 wall | -f 0 RSS | -f 0 wall |
|-----------------------------|-------------|--------------|----------|-----------|
| stock-Java  -Xmx90g         | 12.62 GB    | 36 s         | 15.06 GB | 51 s      |
| stock-Java  -Xmx8g  G1      |  3.81 GB    | 32 s         |  3.65 GB | 43 s      |
| opt-Java    -Xmx8g  G1      |  3.25 GB    | 39 s         |  3.58 GB | 44 s      |
| opt-Java    -Xmx8g  G1 --chunk 50k | 2.04 GB | 39 s     |  2.15 GB | 45 s      |
| **vardictcpp**              | **0.30 GB** | **17 s**     | **0.31 GB** | **19 s** |
| **vardictcpp --chunk 50k**  | **0.02 GB** | **17 s**     | **0.02 GB** | **14 s** |

## Takeaways

- **vardictcpp vs stock-Java (-Xmx90g), -f 0:** 15.06 GB -> 0.31 GB (**49x less memory**),
  51 s -> 19 s (**2.7x faster**). With `--chunk`, 0.02 GB (memory ~independent of region size).
- **Java, byte-identical:** G1 + bounded heap alone cuts stock 15 GB -> 3.6 GB; adding `--chunk`
  reaches 2.15 GB with identical variant content (only the Seg column changes, as with a pre-split BED).
- The `-f 0` regression over `-f 0.01` is modest once the collector is bounded (Java) or absent (C++):
  the dominant term is per-region live set, which `--chunk` and native memory management both bound.

## Count-agreement validation (vardictcpp vs VarDictJava 1.8.3 pileup)

Matched SNV alleles on the same region: **Depth exact 97.7%** (mean abs diff 0.03 reads),
**AltDepth exact 93.3%** (0.14 reads). Residual ±1-read differences are VarDict's paired-read
mate-overlap de-counting, not yet ported (see README "Parity status").
