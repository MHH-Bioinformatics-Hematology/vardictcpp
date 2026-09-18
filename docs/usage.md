# Usage

vardictcpp reads an indexed BAM/CRAM and a reference FASTA, calls variants over a set of regions, and
writes VarDict simple-mode TSV (or VCF with `--vcf`). It accepts VarDict-Java's option syntax and
defaults, so a VarDict-Java command line works unchanged.

The reference FASTA is indexed automatically with htslib `faidx` on first use, so a separate
`samtools faidx` step is not required.

## Quick start

```bash
# a single region
vardictcpp -G ref.fa -b sample.bam -N sample -R chr1:1647000-1648500 -f 0.01

# a BED of target regions (recommended: keeps memory low)
vardictcpp -G ref.fa -b sample.bam -N sample -c 1 -S 2 -E 3 -g 4 -z targets.bed -f 0.01

# eight threads
vardictcpp -G ref.fa -b sample.bam -N sample -c 1 -S 2 -E 3 -g 4 -z targets.bed --th 8
```

## Choosing the regions to call

vardictcpp holds all per-position evidence for the region it is currently processing in memory at
once. Restricting calling to your actual targets keeps peak memory low.

- **A BED of targets (recommended).** Pass the BED as a positional argument and give the column
  indices with `-c` (chromosome), `-S` (start), `-E` (end), and `-g` (gene/name). Add `-z` if the
  BED is 0-based (standard BED).
- **A single region.** `-R chr:start-end` calls one interval directly, without a BED.

!!! warning "Whole chromosomes"
    Passing a whole chromosome to `-R` processes it as one unit and holds all of its per-position data
    in memory at once, which can use a very large amount of RAM. For anything larger than a locus,
    prefer a BED, or split long intervals with `--chunk N` (see [Options](options.md)), which bounds
    peak memory independently of interval length.

## Single-sample calling

```bash
vardictcpp -G ref.fa -b sample.bam -N sample -c 1 -S 2 -E 3 -g 4 -z targets.bed -f 0.01
```

For measurable residual disease (MRD) calling at ultra-low variant allele frequency, set the frequency
threshold to zero so the faintest variant-supporting reads are retained:

```bash
vardictcpp -G ref.fa -b sample.bam -N sample -c 1 -S 2 -E 3 -g 4 -z targets.bed -f 0
```

## Paired (tumor/normal) calling

Provide the two BAMs and the two sample names joined by a pipe. vardictcpp runs the pipeline on both
BAMs and compares them, emitting the somatic STATUS field.

```bash
vardictcpp -G ref.fa -b 'tumor.bam|normal.bam' -N 'tumor|normal' \
  -c 1 -S 2 -E 3 -g 4 -z targets.bed -f 0.01
```

## Multi-threading

`--th N` (alias `--threads`) distributes the target regions across `N` worker threads. Each region is
processed independently and the per-region results are re-assembled into VarDict's deterministic
output order, so the output is byte-identical regardless of the thread count. Threads change only
runtime and peak memory, not the variant calls.

## Output

The default output is VarDict simple-mode TSV, with the same 36-column order as VarDict-Java's
`SimpleOutputVariant`. Passing `--vcf` writes VCF directly (the strand-bias Fisher test and VCF
formatting are done in C++), with one sample column in single-sample mode and tumor and normal columns
plus a somatic STATUS in paired mode.
