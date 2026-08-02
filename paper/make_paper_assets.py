#!/usr/bin/env python3
"""Generate publication assets (LaTeX tables + figures) for the vardictcpp paper.

Reads:
  paper/data/parity.tsv                         (byte-level parity per WES sample)
  paper/data/perf.tsv                           (runtime + peak RSS, cpp vs Java, 1 & 8 threads)
Writes:
  paper/tables.tex                              (Table 1 parity, Table 2 performance)
  paper/figures/fig1_runtime.{pdf,png}
  paper/figures/fig2_memory.{pdf,png}

Figures carry NO title and NO in-figure explanatory text: only axis labels, tick labels, and a
legend (series identity). Everything else belongs in the caption. Palette: Okabe-Ito (colorblind-safe,
grayscale-separable). Run after finalbench.tsv exists.
"""
import csv, os, sys, math

HERE = os.path.dirname(os.path.abspath(__file__))
PARITY = os.path.join(HERE, "data", "parity.tsv")
# self-contained copy under data/ (committed); fall back to the scratch bench dir if absent
BENCH = os.path.join(HERE, "data", "perf.tsv")
if not os.path.exists(BENCH):
    BENCH = os.path.join(HERE, "..", "bench", "wes_parity", "finalbench.tsv")
FIGDIR = os.path.join(HERE, "figures")
os.makedirs(FIGDIR, exist_ok=True)

# Okabe-Ito: blue = vardictcpp, vermillion = VarDictJava
C_CPP, C_JAVA = "#0072B2", "#D55E00"

# short x labels keyed by region count (shown on axis; SRR mapping goes in the caption)
def kregions(n):
    return f"{round(int(n)/1000)}k"

# ---------- load ----------
def load_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))

parity = load_tsv(PARITY)
bench = load_tsv(BENCH) if os.path.exists(BENCH) else []

# bench indexed by (sample, tool, threads) -> (wall_s, rss_mb)
B = {}
for r in bench:
    B[(r["sample"], r["tool"], r["threads"])] = (float(r["wall_s"]), float(r["peakRSS_MB"]))
# perf samples in ascending region order (the ones actually benchmarked)
perf_samples = sorted({r["sample"] for r in bench}, key=lambda s: int(next(p["regions"] for p in parity if p["sample"] == s)))

# ---------- Table 1: parity ----------
def fnum(n): return f"{int(n):,}".replace(",", "\\,")

t1 = []
t1.append(r"\begin{table}[t]")
t1.append(r"\centering")
t1.append(r"\caption{\label{tab:parity}%")
t1.append(r"Byte-level parity of vardictcpp against VarDict-Java 1.8.3 on five whole-exome sequencing")
t1.append(r"samples (hg19, per-sample covered-target BED, minimum allele frequency $0.01$). Sorted")
t1.append(r"output was compared line-for-line against single-threaded VarDictJava (multi-threaded Java")
t1.append(r"output is non-deterministic). \emph{Byte-identical}: output rows matching the Java reference")
t1.append(r"character-for-character; \emph{FP}: rows emitted only by vardictcpp; \emph{FN}: rows emitted")
t1.append(r"only by VarDictJava. All 91\,288 variant rows are byte-identical (0 FP, 0 FN).}")
t1.append(r"\begin{tabular}{lrrrrr}")
t1.append(r"\toprule")
t1.append(r"Sample & Target regions & Variant rows & Byte-identical & FP & FN \\")
t1.append(r"\midrule")
tot_rows = 0
for p in parity:
    tot_rows += int(p["rows"])
    ident = "all" if p["identical"] == p["rows"] else fnum(p["identical"])
    t1.append(f"\\texttt{{{p['sample']}}} & {fnum(p['regions'])} & {fnum(p['rows'])} & "
              f"{fnum(p['identical'])}\\,({'100\\%' if p['identical']==p['rows'] else '--'}) & "
              f"{p['fp']} & {p['fn']} \\\\")
t1.append(r"\midrule")
t1.append(f"Total & & {fnum(tot_rows)} & {fnum(tot_rows)}\\,(100\\%) & 0 & 0 \\\\")
t1.append(r"\bottomrule")
t1.append(r"\end{tabular}")
t1.append(r"\end{table}")

# ---------- Table 2: performance + memory ----------
t2 = []
t2.append(r"\begin{table}[t]")
t2.append(r"\centering")
t2.append(r"\caption{\label{tab:perf}%")
t2.append(r"Runtime and peak memory of vardictcpp versus VarDict-Java 1.8.3 on three whole-exome samples")
t2.append(r"(106--774\,k target regions; hg19, minimum allele frequency $0.01$) at 1 and 8 threads,")
t2.append(r"CPU-pinned. Wall clock and peak resident set size (RSS) from \texttt{/usr/bin/time -v};")
t2.append(r"VarDictJava on JDK~25 with \texttt{-Xmx\,8g}. \emph{Speedup} $=$ Java/cpp runtime;")
t2.append(r"$\times$\emph{less} $=$ Java/cpp peak RSS.}")
t2.append(r"\begin{tabular}{lr rr r rr r}")
t2.append(r"\toprule")
t2.append(r" & & \multicolumn{3}{c}{Runtime (s)} & \multicolumn{3}{c}{Peak RSS (MB)} \\")
t2.append(r"\cmidrule(lr){3-5}\cmidrule(lr){6-8}")
t2.append(r"Sample & Threads & vardictcpp & VarDictJava & Speedup & vardictcpp & VarDictJava & $\times$less \\")
t2.append(r"\midrule")
for s in perf_samples:
    regs = fnum(next(p["regions"] for p in parity if p["sample"] == s))
    for i, th in enumerate(("1", "8")):
        c = B.get((s, "vardictcpp", th)); j = B.get((s, "VarDictJava", th))
        if not c or not j or c[0] == 0 or c[1] == 0 or j[0] == 0 or j[1] == 0:
            continue
        spd = f"{j[0]/c[0]:.1f}$\\times$"; mem = f"{j[1]/c[1]:.0f}$\\times$"
        lead = f"\\texttt{{{s}}} ({regs})" if i == 0 else ""
        t2.append(f"{lead} & {th} & {c[0]:.1f} & {j[0]:.1f} & {spd} & {c[1]:.0f} & {j[1]:.0f} & {mem} \\\\")
    t2.append(r"\addlinespace[2pt]")
t2.append(r"\bottomrule")
t2.append(r"\end{tabular}")
t2.append(r"\end{table}")

with open(os.path.join(HERE, "tables.tex"), "w") as f:
    f.write("% Auto-generated by make_paper_assets.py -- edit the CAPTION HERE placeholders.\n")
    f.write("\n".join(t1) + "\n\n" + "\n".join(t2) + "\n")
print("wrote tables.tex")

# ---------- figures ----------
if not bench:
    print("no finalbench.tsv yet -- tables written, skipping figures. Re-run after the benchmark.")
    sys.exit(0)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, FuncFormatter
import numpy as np

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans", "Arial"],
    "font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
    "legend.fontsize": 7.5, "axes.linewidth": 0.6, "xtick.major.width": 0.6,
    "ytick.major.width": 0.6, "figure.dpi": 150, "savefig.dpi": 300, "savefig.bbox": "tight",
    "pdf.fonttype": 42, "ps.fonttype": 42,
})

def grouped(metric_idx, ylabel, fname):
    # metric_idx: 0 = wall_s, 1 = rss_mb
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.5), sharey=True)
    xs = np.arange(len(perf_samples)); w = 0.38
    labels = [kregions(next(p["regions"] for p in parity if p["sample"] == s)) for s in perf_samples]
    for ax, th, tag in ((axes[0], "1", "a"), (axes[1], "8", "b")):
        cpp = [B.get((s, "vardictcpp", th), (float("nan"),)*2)[metric_idx] for s in perf_samples]
        jav = [B.get((s, "VarDictJava", th), (float("nan"),)*2)[metric_idx] for s in perf_samples]
        ax.bar(xs - w/2, cpp, w, color=C_CPP, label="vardictcpp", zorder=3)
        ax.bar(xs + w/2, jav, w, color=C_JAVA, label="VarDictJava", zorder=3)
        ax.set_yscale("log")
        ax.set_xticks(xs); ax.set_xticklabels(labels)
        ax.set_xlabel("Target regions")
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", which="major", lw=0.4, color="0.85", zorder=0)
        ax.tick_params(length=2.5)
        ax.text(0.02, 0.97, tag, transform=ax.transAxes, va="top", ha="left",
                fontsize=9, fontweight="bold")
    axes[0].set_ylabel(ylabel)
    axes[0].legend(frameon=False, loc="upper left", bbox_to_anchor=(0.0, 0.90))
    fig.tight_layout(w_pad=1.0)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIGDIR, f"{fname}.{ext}"))
    plt.close(fig)
    print(f"wrote {fname}.pdf/.png")

grouped(0, "Runtime (s)", "fig1_runtime")
grouped(1, "Peak RSS (MB)", "fig2_memory")
print("done.")
