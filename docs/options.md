# Options

vardictcpp accepts VarDict-Java 1.8.3's complete option syntax (commons-cli style, including
single-dash multi-character options such as `-th 8`). Every option keeps VarDict's meaning and default,
so existing command lines and wrapper scripts work unchanged. The most commonly used flags are listed
below, with their defaults in parentheses.

## Required inputs

| Flag | Argument | Meaning |
|---|---|---|
| `-G` | FILE | Reference FASTA. Auto-indexed with htslib `faidx`, so no `samtools faidx` step is needed. |
| `-b` | FILE | Input BAM/CRAM. Paired (somatic) mode: `-b 'tumor.bam\|normal.bam'`. |
| `-N` | STR | Sample name. Paired mode: `-N 'tumor\|normal'`. |

## Region selection

| Flag | Argument | Meaning |
|---|---|---|
| `-R` | chr:start-end | Call a single region. Held in memory as one unit, so use a BED or `--chunk` for large intervals. |
| _(BED)_ | path | Positional BED of target intervals to call. |
| `-c` `-S` `-E` `-g` | INT | 1-based BED column indices for chromosome, start, end and gene/name (for example `-c 1 -S 2 -E 3 -g 4`). |
| `-z` | | Treat BED start/end as 0-based (standard BED). |
| `-x` | INT | Extend every region by INT bp up- and downstream, for example `1000` for plus/minus 1 kb (0). |
| `--chunk` | INT | Split any region longer than INT bp into consecutive windows, bounding peak memory. |

## Calling parameters

| Flag | Argument | Meaning |
|---|---|---|
| `-f` | FLOAT | Minimum variant allele frequency (0.01). Set `0` for MRD / ultra-low-VAF calling. |
| `-r` | INT | Minimum alt reads to call a variant (2). |
| `-q` | FLOAT | Base-quality boundary for good/bad base counting (22.5). |
| `-B` | INT | Minimum reads per strand for the strand-bias flag (2). |
| `-O` | FLOAT | Minimum mean mapping quality (0). |
| `-P` | INT | Minimum mean read position for a variant (5). |
| `-o` | FLOAT | Minimum high/low base-quality ratio (1.5). |
| `-X` | INT | Bases inspected past an indel for mismatches (3). |
| `-m` | INT | Skip a read whose mismatches, excluding indel length, exceed this (8). |
| `-I` | INT | Indel-size / large-indel breakpoint search window (50). |
| `-V` | FLOAT | Somatic low-frequency threshold for LOH/somatic gating (0.05, paired mode). |
| `-mfreq` | FLOAT | Variant-frequency threshold for monomer microsatellite regions (0.25). |
| `-nmfreq` | FLOAT | Variant-frequency threshold for non-monomer microsatellite regions (0.1). |
| `-k` | 0\|1 | Local realignment (on). |
| `-t` | | Remove duplicate reads before calling. |

## Output and performance

| Flag | Argument | Meaning |
|---|---|---|
| `--vcf` | | Write VCF directly (in-C++ strand-bias Fisher test and var2vcf). Default output is VarDict simple-mode TSV. |
| `-th` / `--threads` | INT | Region-parallel worker threads (1). Output is byte-identical regardless of thread count. |

!!! note "Full option set"
    vardictcpp parses VarDict-Java's entire option set (62 options) for compatibility. Options that
    drive the pipeline are acted on; the remainder are accepted so that existing VarDict-Java command
    lines are not rejected.
