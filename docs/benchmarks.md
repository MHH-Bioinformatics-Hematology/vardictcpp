# Benchmarks

vardictcpp was benchmarked against VarDict-Java 1.8.3 on five public whole-exome sequencing runs
aligned to hg19, each with its own covered-target BED (11k to 774k regions), at a minimum allele
frequency of 0.01. Wall clock and peak resident set size (RSS) were captured with `/usr/bin/time -v`
on the same host.

## Runtime and peak memory

Single-thread and eight-thread runtimes and peak memory, VarDict-Java on JDK 25 with `-Xmx 8g`
(median of five replicate runs):

| Sample (regions) | Threads | vardictcpp | VarDict-Java | Speedup | vardictcpp RSS | VarDict-Java RSS |
|---|---|---|---|---|---|---|
| SRR15006376 (106k) | 1 | 44.5 s | 181.1 s | 4.1x | 74 MB | 1274 MB |
| SRR15006376 (106k) | 8 | 7.8 s | 85.0 s | 11.0x | 358 MB | 2633 MB |
| SRR15006540 (135k) | 1 | 61.7 s | 226.3 s | 3.7x | 76 MB | 1084 MB |
| SRR15006540 (135k) | 8 | 10.4 s | 96.9 s | 9.3x | 312 MB | 2569 MB |
| SRR8657348 (774k) | 1 | 98.0 s | 473.5 s | 4.8x | 122 MB | 977 MB |
| SRR8657348 (774k) | 8 | 14.8 s | 142.4 s | 9.6x | 252 MB | 2084 MB |

Across the datasets, vardictcpp is 3.7 to 4.8 times faster single-core and 9.3 to 11 times faster on
eight cores, at 7 to 17 times less peak memory.

## Whole-chromosome memory

Memory is highest when the region parameter is set to a whole chromosome, because the per-position
maps then span the entire chromosome. Processing chromosome 1 of SRR8657348 as one work unit,
VarDict-Java peaks at 107.6 GB and 6:10, while vardictcpp holds 10.6 GB and completes in 2:03. Passing
the region as a BED, or splitting it with `--chunk`, keeps memory far lower.

## Accuracy

On all five whole-exome samples, vardictcpp and VarDict-Java produce byte-identical output (the same
variant rows, character for character). On a germline benchmark (chromosome 20 of the Genome in a
Bottle sample HG002, scored against the NIST GIAB v4.2.1 truth set with `rtg vcfeval`), the two
VarDict implementations return identical F1 scores (0.985 for SNVs, 0.801 for indels).

See the manuscript and the Zenodo archive for the full benchmark commands, target BEDs, and
per-replicate measurements.
