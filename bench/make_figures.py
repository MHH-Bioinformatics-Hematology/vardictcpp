#!/usr/bin/env python3
# Publication figures: speed + memory, stock (unfixed) VarDictJava 1.8.3 vs vardictcpp, final binary,
# 3 real WES samples (hg19, bwa mem, covered-target BED, -f 0.01), single- and multi-threaded.
# No titles / no baked-in captions. Vector (PDF/SVG) + 600-dpi PNG.
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, FuncFormatter
import os

OUT = os.path.join(os.path.dirname(__file__), "figures"); os.makedirs(OUT, exist_ok=True)

samples = ["SRR15006540", "SRR15006376", "SRR8657348"]
regions = ["135k targets", "106k targets", "774k targets"]
xlabels = [f"{s}\n{r}" for s, r in zip(samples, regions)]

# wall-clock seconds
java1 = [241.94, 181.88, 502.64];  cpp1 = [58.32, 43.27, 103.17]
java4 = [118.07,  99.99, 191.28];  cpp4 = [15.61, 11.95,  27.67]
# peak RSS (GB)
mjava1 = [1.886, 5.012, 4.171];    mcpp1 = [0.062, 0.073, 0.119]
mjava4 = [2.768, 4.153, 5.396];    mcpp4 = [0.191, 0.198, 0.175]

# ---- shared aesthetics -------------------------------------------------------
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],
    "font.size": 11,
    "axes.linewidth": 0.9,
    "axes.edgecolor": "#3a3f47",
    "xtick.color": "#3a3f47", "ytick.color": "#3a3f47",
    "axes.labelcolor": "#22262b",
    "svg.fonttype": "none",
})
C_JAVA_1, C_JAVA_4 = "#B9C4CE", "#5B6B7B"   # baseline: muted cool grey-blue (light=1t, dark=4t)
C_CPP_1,  C_CPP_4  = "#7FD8C4", "#0E9C7E"   # ours: vivid teal/green
EDGE = "#2b2f36"

def grouped(ax, series, colors, labels, value_fmt, ymin, ymax):
    x = np.arange(len(samples)); w = 0.20; offs = [-1.5*w, -0.5*w, 0.5*w, 1.5*w]
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, which="major", color="#dfe3e8", lw=0.8)
    ax.yaxis.grid(True, which="minor", color="#eef1f4", lw=0.5)
    bars = []
    for vals, c, off in zip(series, colors, offs):
        b = ax.bar(x + off, vals, w, color=c, edgecolor=EDGE, linewidth=0.6, zorder=3)
        bars.append(b)
        for xi, v in zip(x + off, vals):
            ax.annotate(value_fmt(v), (xi, v), textcoords="offset points", xytext=(0, 2.5),
                        ha="center", va="bottom", fontsize=7.3, color="#22262b", rotation=90)
    ax.set_yscale("log")
    ax.set_ylim(ymin, ymax)
    ax.set_xticks(x); ax.set_xticklabels(xlabels, fontsize=10)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    ax.tick_params(axis="x", length=0)
    ax.legend(bars, labels, frameon=False, fontsize=9.2, ncol=2,
              loc="upper center", bbox_to_anchor=(0.5, 1.14), columnspacing=1.6, handlelength=1.3)
    return x

def save(fig, name):
    for ext, kw in (("pdf", {}), ("svg", {}), ("png", {"dpi": 600})):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"), bbox_inches="tight", **kw)

# ---- speed -------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.4, 4.3))
grouped(ax, [java1, cpp1, java4, cpp4], [C_JAVA_1, C_CPP_1, C_JAVA_4, C_CPP_4],
        ["VarDictJava  (1 thread)", "vardictcpp  (1 thread)",
         "VarDictJava  (4 threads)", "vardictcpp  (4 threads)"],
        lambda v: f"{v:.0f}s" if v >= 10 else f"{v:.1f}s", 6, 1500)
ax.set_ylabel("Wall-clock time  (s, log scale)", fontsize=11.5)
ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{y:g}"))
fig.tight_layout(); save(fig, "fig_speed"); plt.close(fig)

# ---- memory ------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.4, 4.3))
grouped(ax, [mjava1, mcpp1, mjava4, mcpp4], [C_JAVA_1, C_CPP_1, C_JAVA_4, C_CPP_4],
        ["VarDictJava  (1 thread)", "vardictcpp  (1 thread)",
         "VarDictJava  (4 threads)", "vardictcpp  (4 threads)"],
        lambda v: f"{v:.2f}", 0.03, 22)
ax.set_ylabel("Peak resident memory  (GB, log scale)", fontsize=11.5)
ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{y:g}"))
fig.tight_layout(); save(fig, "fig_memory"); plt.close(fig)

# ---- stats to stdout ---------------------------------------------------------
def row(n, j, c): return f"{n:14s} Java {j:8.1f}  C++ {c:7.1f}  ->  {j/c:5.2f}x"
print("SPEED (wall s), speedup = Java/C++")
for i, s in enumerate(samples):
    print(" ", row(f"{s} 1t", java1[i], cpp1[i])); print(" ", row(f"{s} 4t", java4[i], cpp4[i]))
print(f"  mean speedup 1t: {np.mean([j/c for j,c in zip(java1,cpp1)]):.2f}x   "
      f"4t: {np.mean([j/c for j,c in zip(java4,cpp4)]):.2f}x")
print("\nMEMORY (peak GB), reduction = Java/C++")
for i, s in enumerate(samples):
    print(f"  {s} 1t: Java {mjava1[i]:.3f}  C++ {mcpp1[i]:.3f} -> {mjava1[i]/mcpp1[i]:5.1f}x   "
          f"4t: Java {mjava4[i]:.3f}  C++ {mcpp4[i]:.3f} -> {mjava4[i]/mcpp4[i]:5.1f}x")
print(f"  memory reduction range: {min(mjava4[i]/mcpp4[i] for i in range(3)):.0f}-"
      f"{max(mjava1[i]/mcpp1[i] for i in range(3)):.0f}x")
print(f"\nfigures written to {OUT}/  (pdf, svg, png)")
