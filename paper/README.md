# Publication assets — vardictcpp (GigaByte software paper)

Self-contained source for the manuscript figures and tables. Regenerate everything with:

```bash
python3 make_paper_assets.py    # needs matplotlib + numpy
```

## Contents

| File | What it is |
|---|---|
| `outline.md` | Section-by-section manuscript outline (bullet points), GigaByte Technical Release. `[VERIFY]`/`[TODO]` marks items to confirm before submission. |
| `data/parity.tsv` | Byte-level parity per WES sample (regions, rows, byte-identical, FP, FN). Source for Table 1. |
| `data/perf.tsv` | Runtime + peak RSS, vardictcpp vs VarDictJava, 1 & 8 threads (clean CPU-pinned runs). Source for Table 2 and the performance figures. Add a `repbench.tsv` (with a `rep` column) to draw replicate box plots. |
| `data/accuracy.tsv` | GIAB HG002 chr20 exome precision/recall/F1 by class (All/SNV/Indel). Source for Table 3. |
| `data/per_region.tsv` | Per-1Mb-window precision/recall/F1 (both implementations). Source for the accuracy box plots. |
| `data/benchmark.sh` | Exact performance-benchmark script (CPU-pinned, `/usr/bin/time -v`, JDK 25, `-Xmx 8g`). |
| `data/giab_accuracy.sh` | Exact accuracy pipeline (chr20 fetch, calling, `var2vcf`, `rtg vcfeval`). |
| `data/strandbias.py` | Pure-Python reproduction of VarDict's `teststrandbias.R` (no R/scipy needed). |
| `data/per_region_metrics.py` | Computes per-window precision/recall/F1 from `rtg vcfeval` output. |
| `make_paper_assets.py` | Generates `tables.tex` + `figures/fig1_runtime` + `fig2_memory` + `fig3_accuracy`. |
| `tables.tex` | Table 1 (parity), 2 (performance), 3 (accuracy), captions included. |
| `figures.tex` | `figure` blocks with full captions (figures carry no in-figure prose). |
| `figures/` | Publication-ready PDFs (vector) + PNG previews. |

## Figures (all box plots, no bars)

- **Fig 1 / Fig 2** — runtime and peak-RSS, vardictcpp vs VarDictJava, 1- and 8-thread panels, log-y; each clean CPU-pinned run is a marker (boxes appear where ≥3 replicates are supplied).
- **Fig 3** — per-region precision/recall/F1 on the GIAB HG002 chr20 exome, SNV vs indel; the vardictcpp and VarDictJava boxes coincide exactly (byte-identical VCFs).

## Conventions

- No titles, no in-figure explanatory text: only axis labels, tick labels, a legend, and small `a`/`b`/`c` panel tags. Everything descriptive lives in the caption (`figures.tex`).
- Log-scaled y-axis for runtime and memory; boxes = median + IQR + 1.5×IQR whiskers, individual runs/windows overplotted.
- Okabe-Ito colorblind-safe palette: blue `#0072B2` = vardictcpp, vermillion `#D55E00` = VarDictJava.
- Vector PDF at 300 dpi, Type-42 (editable) fonts.
