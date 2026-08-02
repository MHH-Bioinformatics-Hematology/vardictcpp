# Publication assets — vardictcpp (GigaByte software paper)

Self-contained source for the manuscript figures and tables. Regenerate everything with:

```bash
python3 make_paper_assets.py    # needs matplotlib + numpy
```

## Contents

| File | What it is |
|---|---|
| `outline.md` | Section-by-section manuscript outline (bullet points), GigaByte Technical Release. `[VERIFY]` marks items to confirm before submission. |
| `data/parity.tsv` | Byte-level parity per WES sample (regions, rows, byte-identical, FP, FN). Source for Table 1. |
| `data/perf.tsv` | Runtime + peak RSS, vardictcpp vs VarDictJava 1.8.3, 1 and 8 threads. Source for Table 2 and both figures. |
| `data/benchmark.sh` | The exact benchmark script that produced `perf.tsv` (CPU-pinned, `/usr/bin/time -v`, JDK 25, `-Xmx 8g`). |
| `make_paper_assets.py` | Generates `tables.tex` + `figures/fig1_runtime.{pdf,png}` + `figures/fig2_memory.{pdf,png}`. |
| `tables.tex` | Table 1 (`tab:parity`) and Table 2 (`tab:perf`), captions included. |
| `figures.tex` | `figure` blocks with full captions (the figures carry no in-figure prose — all text is in the captions). |
| `figures/` | Publication-ready PDFs (vector) + PNG previews. |

## Figure conventions

- No titles, no in-figure explanatory text: only axis labels, tick labels, a legend, and small `a`/`b` panel tags. Everything descriptive lives in the caption (`figures.tex`).
- Log-scaled y-axis (runtime and memory span an order of magnitude).
- Okabe-Ito colorblind-safe palette: blue `#0072B2` = vardictcpp, vermillion `#D55E00` = VarDictJava.
- Vector PDF at 300 dpi, Type-42 (editable) fonts.
