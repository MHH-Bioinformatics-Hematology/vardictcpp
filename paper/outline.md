# vardictcpp — manuscript outline (GigaByte / GigaScience software paper)

Target: **GigaByte** (GigaScience Press) *Technical Release* — a software tool paper. Keep it concise,
reproducibility-first; GigaByte expects a runnable artifact, an archived release, and a data/parity
statement. Fill every `[VERIFY]` before submission.

---

## Title (working)
- "vardictcpp: a byte-identical, memory-lean C++ reimplementation of the VarDict variant caller"
- Alternatives emphasising the contribution: "…a drop-in C++ port of VarDict with verified byte-level
  parity and an order-of-magnitude lower memory footprint."

## Abstract (structured, ~150–200 words)
- **Background**: VarDict (Perl) / VarDictJava is a widely used amplicon- and WES/WGS-capable SNV+indel
  caller; the JVM implementation is memory-heavy and its startup/GC dominate small jobs.
- **Findings**: a from-scratch C++17 reimplementation (htslib) that reproduces VarDictJava 1.8.3 output
  **byte-for-byte** across all modes and on five whole-exome samples (11k–774k target regions), while
  running **~3.7–4.8× faster single-core / ~9–11× at 8 threads** and using **~8–17× less peak RAM**
  (74–122 MB vs ~1.0–1.3 GB single-threaded; ~2.1–2.6 GB at 8 threads). Against an independent
  ground-truth set (GIAB HG002 exome, NIST v4.2.1) the two implementations are byte-identical, so they
  share exactly the same accuracy (SNV F$_1$ 0.985, indel F$_1$ 0.799 in confident regions).
- **Conclusion**: a verified drop-in replacement enabling VarDict on memory-constrained and
  high-throughput settings; every divergence encountered during development was a port bug corrected
  toward the reference (no upstream-bug divergences), and the two latent bugs found in the process are
  reported.

## Keywords
- variant calling; VarDict; SNV; indel; structural variants; C++; htslib; reproducibility; byte-level
  parity; memory footprint

---

## Background / Context
- Variant calling in clinical/MRD and cancer pipelines needs speed, low memory, and *reproducibility*.
- VarDict's algorithm (per-position pileup + local realignment of soft-clipped reads for indels/SV,
  amplicon-bias-aware, somatic tumor–normal) — cite Lai et al. 2016 (NAR, doi:10.1093/nar/gkw227).
- The Perl and Java implementations: upstream declares 1.8.3 the final VarDictJava version and no longer
  actively maintains it — motivates an independent, maintained, resource-lean implementation.
- Problem statement: (i) JVM memory/GC overhead (up to tens of GB, GC hoards toward `-Xmx`); (ii) JVM
  startup cost on many small regions; (iii) reproducibility risk when reimplementing — most rewrites are
  "concordant" (≈ set-equal), not byte-identical, so downstream diffs are hard to attribute.
- **Contribution**: a C++ port validated to *byte-level* parity (not just variant-set concordance), with
  a resource benchmark; the byte-level bar is what makes it a trustworthy drop-in.

## Methods / Implementation
- **Architecture**: mirror the VarDictJava package layout module-for-module so parity is auditable
  (config, region, reference, cigar parse, cigar modifier, realignment, structural variants, to-vars /
  post-process, printers, somatic, amplicon, fisher). One C++ file per Java module; table of the mapping.
- **I/O**: htslib for BAM/CRAM (`sam_read1`/`bam1_t`, CIGAR ops) and FASTA (`faidx`); CLI accepts
  VarDictJava's full commons-cli option set (single-dash multi-char options, `-th`, `-VS`, `-mfreq`…).
- **Numeric parity landmines** (call these out — they are the interesting engineering): Perl-semantics
  `substr`/`charAt` negative indices; `DecimalFormat`/round-half-to-even; the reference "genuine window"
  map (disjoint, `SEED_1`-truncated) vs a naive contiguous buffer; the CigarModifier ordering; somatic
  `combineAnalysis` re-entrancy on a merged tumor+normal region.
- **Modes ported**: simple, amplicon, `--fisher`, structural variants (split-read + pair-assisted
  `<INV>`; split-read + discordant `<DEL>`; tandem-duplication `<DUP>` incl. inter-chromosomal fusion
  clusters + the SOFTP2SV guard), paired somatic incl. `combineAnalysis`, splice-junction handling.
- **Memory design**: per-region streaming (one record at a time, no read list); reference stored as
  *disjoint* windows matching Java's reference map (a far realignment breakpoint never gap-fills a
  multi-Mbp span); an O(n) rolling k-mer seed index (vs O(n·k)); portable SIMD (SSE2/NEON/scalar) for the
  base transforms.
- **Verification protocol** (the methodological core): define parity as `diff == 0` on the sorted TSV vs
  single-threaded VarDictJava 1.8.3 (note: multi-threaded Java output is non-deterministic — always
  compare against `-th 1`). Curated golden fixtures per mode (CI-gated) + five public WES samples.
  Every residual row triaged read-by-read (`samtools view` + VarDict `-y` trace) against **both** Java
  and the original Perl to classify port-bug vs upstream-bug; port bugs fixed, upstream bugs would be
  documented as intentional divergences.
- **Ground-truth accuracy** (independent of the Java reference): call VarDict on a benchmark sample
  with a public high-confidence truth set (GIAB HG002, Agilent SureSelect v5 exome, chromosome 20) and
  score with `rtg vcfeval` over callable confident regions (exome depth ≥ 20; germline operating point
  = PASS and AF ≥ 0.2). Report per-region precision/recall/F1 distributions. Because vardictcpp and
  VarDictJava emit **byte-identical VCFs on this sample too**, this measures VarDict's accuracy while
  showing the port inherits it exactly. hg19 primary-assembly sequence equals GRCh37 for chr1–22,X, so
  the local hg19 reference (contig renamed to match the truth) is reused without a separate download.

## Results / Validation
- **Byte-level parity** (Table 1): all five WES samples `diff == 0` (0 FP / 0 FN), plus per-mode golden
  fixtures. State the total rows compared (≈ 91k variant rows across the five samples).
- **Correctness triage outcome**: of ~40 residual rows across the five samples, *every one* was a port
  bug corrected toward Java/Perl — `DIVERGENCES.md` has no open entries. This is the headline
  reproducibility result.
- **Two latent bugs surfaced** (worth a short subsection — value beyond the port): (i) an output-buffer
  overflow that silently dropped large-allele variant lines genome-wide [note: whether this is
  vardictcpp-only or also affects a downstream consumer — [VERIFY]]; (ii) a soft-clip read filter that in
  vardictcpp was applied to the pre-modification CIGAR, manufacturing false indels (matched Java by
  moving it post-`modifyCigar`).
- **Ground-truth accuracy** (Table 3, Fig 3): on GIAB HG002 chr20 exome vs NIST v4.2.1, VarDict scores
  SNV precision/recall/F1 = 0.996/0.974/0.985 and indel 0.877/0.734/0.799 in callable confident regions
  (2\,443 truth variants). vardictcpp and VarDictJava emit **byte-identical VCFs** here as well
  (8\,244 variants, 0 differing lines), so the per-region precision/recall/F1 box plots for the two
  implementations coincide exactly — the port preserves accuracy with no deviation.
- **Performance** (Table 2, Fig 1): wall-clock cpp vs Java, single-core and 8-thread, 3 WES samples,
  box plots over 5 clean (contention-filtered) replicate runs.
- **Memory** (Table 2, Fig 2): peak RSS cpp vs Java; ~8–17× reduction single-threaded (74–122 MB vs
  ~1.0–1.3 GB), ~7–9× at 8 threads; note the reference-window design is
  what bounds it (a naive contiguous reference gap-fill inflated peak RSS ~25× — mention as a design
  lesson, Fig 3 optional).
- **Scaling** (Fig 3 optional): cpp thread scaling on the 774k-region sample.

## Availability & requirements (GigaByte requires this block)
- Project name: vardictcpp
- Project home page: `https://github.com/[ORG]/vardictcpp` [VERIFY the canonical URL]
- Archived release: deposit a tagged release + Zenodo DOI [TODO], and/or GigaDB.
- Operating systems: Linux, macOS (x86-64 and ARM/Apple Silicon — CI matrix)
- Programming language: C++17
- Other requirements: CMake ≥ 3.15, htslib
- License: MIT (retaining the upstream AstraZeneca-NGS copyright)
- RRID / bio.tools ID: register [TODO]

## Reuse potential
- Drop-in for VarDict in existing pipelines (same TSV → `teststrandbias.R`/`var2vcf`), so no downstream
  change; the Galaxy wrapper (this project) exposes it.
- Enables VarDict on memory-constrained nodes and large cohorts (the memory + startup wins compound).
- The byte-parity harness is itself reusable as a template for validating other reimplementations.

## Data availability / reproducibility
- Input WES BAMs: SRA run accessions [Table 1]; reference hg19. Provide the exact commands and the
  covered-target BEDs. Archive the five Java `-th 1` golden TSVs + the comparison scripts as a GigaDB
  dataset [TODO] so parity is independently re-checkable.
- **Sample provenance** (confirmed via ENA): SRR15006375/376/386/540 are targeted-capture AML samples
  from BioProject PRJNA742684 (adult AML bone marrow, Illumina MiSeq); SRR8657348 is the MV4-11 AML
  cell line (BioProject PRJNA523380, CCLE targeted capture).
- **Accuracy benchmark data**: GIAB HG002 (NA24385) Agilent SureSelect v5 exome BAM
  (OsloUniversityHospital_Exome, GRCh37) + NIST v4.2.1 benchmark VCF/BED. Provide the chr20 extraction,
  reference-contig renaming, calling, `var2vcf` + strand-bias steps, and `rtg vcfeval` commands
  (`bench/giab_acc/`). Truth and confident BED are public; archive the derived call VCFs + scripts.

## Declarations
- Competing interests; funding; author contributions; acknowledgements (VarDict/VarDictJava authors —
  AstraZeneca-NGS; cite the original paper).

---

## Figure & table plan (files in `paper/`)
- **Table 1** (`tables.tex`, `tab:parity`): per-sample byte-level parity — sample, source, target
  regions, variant rows, byte-identical, FP, FN; plus a mode-fixture parity row block.
- **Table 2** (`tables.tex`, `tab:perf`): runtime + peak RSS (median of 5 clean replicates), cpp vs
  Java, 1 and 8 threads, per sample, with speedup and memory-reduction factors.
- **Table 3** (`tables.tex`, `tab:accuracy`): GIAB HG002 chr20 exome accuracy (All/SNV/Indel): truth,
  TP, FP, FN, precision, recall, F1 (identical for both implementations).
- **Fig 1** (`figures/fig1_runtime.pdf`): box plots (5 clean replicates), wall-clock runtime, cpp vs
  Java, single-core and 8-thread panels, 3 WES samples. Log y. No in-figure text.
- **Fig 2** (`figures/fig2_memory.pdf`): box plots, peak RSS, cpp vs Java, log y. No title.
- **Fig 3** (`figures/fig3_accuracy.pdf`): per-region precision/recall/F1 box plots (SNV vs indel),
  cpp vs Java coinciding exactly. No title.

## To finalise before writing
- [x] Confirm sample provenance (BioProject/GSE; CCLE identity of SRR8657348) — done via ENA:
      PRJNA742684 (AML WES) and PRJNA523380 (MV4-11 / CCLE).
- [ ] Confirm canonical GitHub URL, license header, and cut a tagged release + Zenodo DOI.
- [ ] Decide the scope claim wording: "byte-identical on the five WES samples tested + all mode
      fixtures" (accurate) — do **not** claim "all possible WES data" (an asymptote).
- [ ] Check whether the buffer-overflow bug also affects VarDictJava (would strengthen the bug report).
- [ ] Accuracy caveat wording: the GIAB benchmark is chr20 only (bounds compute); note this and that it
      characterises VarDict, inherited byte-identically by the port — not a novel-accuracy claim.
