# vardictcpp 0.1.0

First release of **vardictcpp** — a C++17 port of [VarDict](https://github.com/AstraZeneca-NGS/VarDictJava)
(AstraZeneca-NGS VarDictJava 1.8.3), built on htslib. It reproduces VarDict's variant calls while running
several times faster and using an order of magnitude less memory.

## Highlights
- **All five modes ported and byte-identical** to VarDictJava 1.8.3 on the test goldens: simple, `--fisher`,
  somatic (paired), amplicon, and structural variants (`<DEL>` discordant-pair + `<INV>` split-read).
- **Faster and far leaner** than stock VarDictJava (see below).
- **Streaming, region-parallel** (`-th`) with ordered output identical to single-threaded; optional
  `--chunk` bounds memory on very large regions.
- Drop-in CLI: accepts VarDict's option syntax; VarDict simple-mode TSV column order.

## Performance vs stock VarDictJava 1.8.3
Real whole-exome data (3 public SRA runs aligned to hg19, covered-target BED, `-f 0.01`), peak RSS via
`/usr/bin/time -v`:

| metric | single thread | 4 threads |
|---|---|---|
| **Speed** (wall, Java ÷ C++) | **4.4× faster** (mean) | **7.6× faster** (mean) |
| per-sample speedup range | 4.2–4.9× | 6.9–8.4× |
| **Peak memory** (Java ÷ C++) | **30–69× less** | **14–31× less** |
| absolute peak RSS | C++ 60–120 MB vs Java 1.9–5.0 GB | C++ ~0.2 GB vs Java 2.8–5.4 GB |

On the CI fixture the same effect is visible in miniature: **~10× faster, ~50× less peak memory**. The
memory advantage grows with region size — on a 1 Mb region at `-f 0` VarDict retains every covered
base × allele, where the native implementation stays flat instead of scaling into multiple GB.

## Correctness
- **Byte-identical** to VarDictJava 1.8.3 across all five modes on the curated goldens (0 FP / 0 FN),
  enforced in CI on every push (gcc + clang).
- On noisy real whole-exome data, **~95% of variant rows are byte-identical**, with a small residual of
  a few dozen discordant calls per sample (tens of FP/FN out of 16k–54k) in the known edge-case classes
  (structural / large-indel representation, distributed indel coverage, homopolymer/MNV). These are not
  yet zero on full WES; the curated-golden parity is exact.

## Build & test
- Requires a C++17 compiler, CMake ≥ 3.15, and htslib:
  `cmake -S . -B build -DHTSLIB_ROOT=$CONDA_PREFIX && cmake --build build -j`
- `bash test/run_tests.sh` — parity against a self-contained VarDictJava golden fixture.
- GitHub Actions CI: gcc/clang build + parity test, plus a runtime & memory benchmark vs VarDictJava.

## Known limitations
- Real-WES parity is ~95% byte-identical (see Correctness); the remaining SV/large-indel/coverage edge
  cases are the roadmap.
- Splicing mode and differing-sample somatic classification (StrongSomatic/LOH) are not yet validated.
