# vardictcpp 1

First release of **vardictcpp** — a C++17 port of [VarDict](https://github.com/AstraZeneca-NGS/VarDictJava)
(AstraZeneca-NGS VarDictJava 1.8.3), built on htslib. It reproduces VarDict's variant calls while running
several times faster and using an order of magnitude less memory.

## Highlights
- **Simple mode is byte-identical** to VarDictJava 1.8.3 on the curated golden (0 FP / 0 FN), enforced
  in CI. **`--fisher`, amplicon, and structural variants** (`<DEL>` discordant-pair, `<INV>` split-read)
  are ported and enabled; **paired somatic** runs the pipeline on both BAMs and compares them (all 56
  rows byte-identical to Java on the test pair, incl. `combineAnalysis`). See *Known limitations*.
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
- **Simple mode is byte-identical** to VarDictJava 1.8.3 on the curated golden (0 FP / 0 FN), enforced
  in CI on every push (gcc + clang; the CI parity test covers simple mode).
- On real whole-exome data, **byte-identical to Java on all five samples tested**, zero differing lines
  vs single-threaded VarDict-Java 1.8.3: SRR15006386 (4014/4014), SRR15006375 (4056/4056), SRR15006376
  (12739/12739), SRR15006540 (16031/16031), and the 774k-region CCLE run SRR8657348 (54448/54448).
  Reaching this drove the full structural-variant subsystems (`<INV>`/`<DEL>`/`<DUP>` incl. discordant
  pairs and inter-chromosomal fusion clusters), large-indel coverage reloads, the reference `SEED_1`
  extent truncation, and the CigarModifier soft-clip/homopolymer fixes. Divergences were triaged for
  *correctness* against both Java and the original Perl: every one proved to be a cpp bug fixed toward
  the reference, so `docs/DIVERGENCES.md` has no open entries. (Other WES samples are not exhaustively
  verified; new data may surface further edge cases.)
- **Performance/memory** on those samples: ~3.5–4.7x faster single-core, ~10–12x at 8 threads, and ~14x
  less peak RAM (77–122 MB vs Java's 1.1–1.7 GB).

## Build & test
- Requires a C++17 compiler, CMake ≥ 3.15, and htslib:
  `cmake -S . -B build -DHTSLIB_ROOT=$CONDA_PREFIX && cmake --build build -j`
- `bash test/run_tests.sh` — parity against a self-contained VarDictJava golden fixture.
- GitHub Actions CI: gcc/clang build + parity test, plus a runtime & memory benchmark vs VarDictJava.

## Known limitations
- **Paired somatic** runs the pipeline on both BAMs and compares them, **byte-identical to Java** on the
  test tumor|normal pair (all 56 rows), including `combineAnalysis` (the merged `bam1+bam2` refinement),
  verified on a fixture that provably triggers it.
- Real-WES parity is byte-identical on all five verified samples (see Correctness); other samples are
  not exhaustively verified.
- Splice junctions are handled (N-op intron spans reject splice-junction deletions in isGoodVar,
  verified against Java on a synthetic spliced fixture).
