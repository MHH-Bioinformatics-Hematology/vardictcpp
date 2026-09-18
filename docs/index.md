# vardictcpp

**vardictcpp** is a standalone C++17/htslib reimplementation of the
[VarDict](https://github.com/AstraZeneca-NGS/VarDictJava) variant caller. It produces output
byte-for-byte identical to VarDict-Java 1.8.3 on the evaluated whole-exome samples, at a fraction of
the runtime and peak memory, in single-sample and paired (tumor/normal) modes.

## Why vardictcpp

- **Same calls.** vardictcpp reproduces VarDict's calling algorithm; on five whole-exome samples the
  output is byte-identical to VarDict-Java 1.8.3.
- **Faster.** About 4x faster single-core and 9-11x faster on eight cores than VarDict-Java.
- **Much smaller memory footprint.** Up to 17x less peak memory, which lets more calling jobs run
  concurrently on the same node. This matters for ultra-low-VAF measurable residual disease (MRD)
  calling, where the allele-frequency threshold is set to zero.
- **Drop-in.** It accepts VarDict-Java's option syntax and default values and writes the same output,
  so existing pipelines and wrapper scripts work unchanged.

## At a glance

```bash
# install
conda install -c bioconda vardictcpp

# single-sample calling over a BED of targets
vardictcpp -G ref.fa -b sample.bam -N sample -c 1 -S 2 -E 3 -g 4 -z targets.bed -f 0.01

# paired tumor/normal (somatic) calling
vardictcpp -G ref.fa -b 'tumor.bam|normal.bam' -N 'tumor|normal' -c 1 -S 2 -E 3 -g 4 -z targets.bed
```

## Where to go next

- [Installation](installation.md): Bioconda, build from source, Galaxy.
- [Usage](usage.md): quick start, region modes, single and paired calling, output.
- [Options](options.md): the full option reference with defaults.
- [Benchmarks](benchmarks.md): runtime and memory versus VarDict-Java.
- [Support and citation](support.md): where to get help and how to cite.
