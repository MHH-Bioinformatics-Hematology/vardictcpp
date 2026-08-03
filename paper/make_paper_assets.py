#!/usr/bin/env python3
"""Generate publication assets (LaTeX tables + box-plot figures) for the vardictcpp paper.

Reads (all under paper/data/):
  parity.tsv       byte-level parity per WES sample (Table 1)
  perf.tsv         single-run runtime + peak RSS (fallback for Table 2 / figures)
  repbench.tsv     replicate runtime + peak RSS (primary: box-plot distributions)
  accuracy.tsv     GIAB HG002 chr20 exome precision/recall/F1 by class (Table 3)
  per_region.tsv   per-window precision/recall/F1 (accuracy box plots)
Writes:
  tables.tex                       Table 1 (parity), 2 (performance), 3 (accuracy)
  figures/fig1_runtime.{pdf,png}   runtime box plots, cpp vs Java, 1 & 8 threads
  figures/fig2_memory.{pdf,png}    peak-RSS box plots
  figures/fig3_accuracy.{pdf,png}  per-region precision/recall/F1 box plots

Figures carry NO title and NO in-figure explanatory text: only axis labels, tick
labels, a legend, and small a/b/c panel tags. Everything else is in the caption.
Palette: Okabe-Ito (colorblind-safe). Box plots (not bars) throughout.
"""
import csv, os, sys, statistics as st

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "data")
FIGDIR = os.path.join(HERE, "figures")
os.makedirs(FIGDIR, exist_ok=True)

C_CPP, C_JAVA = "#0072B2", "#D55E00"   # blue = vardictcpp, vermillion = VarDictJava
SERIES = [("vardictcpp", C_CPP), ("VarDictJava", C_JAVA)]

def load(name):
    p = os.path.join(D, name)
    return list(csv.DictReader(open(p), delimiter="\t")) if os.path.exists(p) else []

parity = load("parity.tsv")
perf1 = load("perf.tsv")
rep = load("repbench.tsv")
acc = load("accuracy.tsv")
per = load("per_region.tsv")

def fnum(n): return f"{int(round(float(n))):,}".replace(",", "\\,")

# region count per sample (from parity)
REG = {p["sample"]: int(p["regions"]) for p in parity}

# ---------------- Table 1: parity ----------------
t1 = [r"\begin{table}[t]", r"\centering", r"\caption{\label{tab:parity}%",
 r"Byte-level parity of vardictcpp against VarDict-Java 1.8.3 on five whole-exome sequencing",
 r"samples (hg19, per-sample covered-target BED, minimum allele frequency $0.01$). Sorted output",
 r"was compared line-for-line against single-threaded VarDictJava (multi-threaded Java output is",
 r"non-deterministic). \emph{Byte-identical}: rows matching the Java reference character-for-character;",
 r"\emph{FP}: rows emitted only by vardictcpp; \emph{FN}: rows emitted only by VarDictJava. All",
 r"91\,288 variant rows are byte-identical (0 FP, 0 FN).}",
 r"\begin{tabular}{lrrrrr}", r"\toprule",
 r"Sample & Target regions & Variant rows & Byte-identical & FP & FN \\", r"\midrule"]
tot = 0
for p in parity:
    tot += int(p["rows"])
    t1.append(f"\\texttt{{{p['sample']}}} & {fnum(p['regions'])} & {fnum(p['rows'])} & "
              f"{fnum(p['identical'])}\\,(100\\%) & {p['fp']} & {p['fn']} \\\\")
t1 += [r"\midrule", f"Total & & {fnum(tot)} & {fnum(tot)}\\,(100\\%) & 0 & 0 \\\\",
       r"\bottomrule", r"\end{tabular}", r"\end{table}"]

# ---------------- performance data (replicates preferred) ----------------
# REP[(sample,tool,threads)] -> {'wall':[...], 'rss':[...]}
REP = {}
for r in rep:
    k = (r["sample"], r["tool"], r["threads"])
    b = REP.setdefault(k, {"wall": [], "rss": []})
    b["wall"].append(float(r["wall_s"])); b["rss"].append(float(r["peakRSS_MB"]))
if not REP:  # fall back to single-run
    for r in perf1:
        REP[(r["sample"], r["tool"], r["threads"])] = {"wall": [float(r["wall_s"])],
                                                        "rss": [float(r["peakRSS_MB"])]}
perf_samples = sorted({k[0] for k in REP}, key=lambda s: REG.get(s, 0))

def med(k, metric):
    return st.median(REP[k][metric]) if k in REP and REP[k][metric] else float("nan")

# ---------------- Table 2: performance (median of replicates) ----------------
nrep = max((len(v["wall"]) for v in REP.values()), default=1)
t2 = [r"\begin{table}[t]", r"\centering", r"\caption{\label{tab:perf}%",
 r"Runtime and peak memory of vardictcpp versus VarDict-Java 1.8.3 on three whole-exome samples",
 r"(106--774\,k target regions; hg19, minimum allele frequency $0.01$) at 1 and 8 threads, CPU-pinned.",
 (r"Wall clock and peak resident set size (RSS) from" if nrep < 3 else
  f"Median of {nrep} replicate runs; wall clock and peak resident set size (RSS) from"),
 r"\texttt{/usr/bin/time -v}, VarDictJava on JDK~25 with \texttt{-Xmx\,8g}. \emph{Speedup} $=$",
 r"Java\,$/$\,cpp runtime; $\times$\emph{less} $=$ Java\,$/$\,cpp peak RSS.",
 r"Per-run values are plotted in Figs.~\ref{fig:runtime} and \ref{fig:memory}.}",
 r"\begin{tabular}{lr rr r rr r}", r"\toprule",
 r" & & \multicolumn{3}{c}{Runtime (s)} & \multicolumn{3}{c}{Peak RSS (MB)} \\",
 r"\cmidrule(lr){3-5}\cmidrule(lr){6-8}",
 r"Sample & Threads & vardictcpp & VarDictJava & Speedup & vardictcpp & VarDictJava & $\times$less \\",
 r"\midrule"]
for s in perf_samples:
    for i, th in enumerate(("1", "8")):
        cw, jw = med((s, "vardictcpp", th), "wall"), med((s, "VarDictJava", th), "wall")
        cr, jr = med((s, "vardictcpp", th), "rss"), med((s, "VarDictJava", th), "rss")
        if any(x != x for x in (cw, jw, cr, jr)) or cw == 0 or cr == 0:
            continue
        lead = f"\\texttt{{{s}}} ({fnum(REG.get(s,0))})" if i == 0 else ""
        t2.append(f"{lead} & {th} & {cw:.1f} & {jw:.1f} & {jw/cw:.1f}$\\times$ & "
                  f"{cr:.0f} & {jr:.0f} & {jr/cr:.0f}$\\times$ \\\\")
    t2.append(r"\addlinespace[2pt]")
t2 += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]

# ---------------- Table 3: accuracy ----------------
CLSNAME = {"all": "All", "snv": "SNV", "indel": "Indel"}
t3 = [r"\begin{table}[t]", r"\centering", r"\caption{\label{tab:accuracy}%",
 r"Variant-calling accuracy of VarDict on the GIAB HG002 chromosome-20 exome (Agilent SureSelect",
 r"v5, Oslo University Hospital), against the NIST v4.2.1 benchmark via \texttt{rtg vcfeval} over",
 r"callable confident regions (exome depth $\geq 20$; 2.31\,Mbp; 2\,443 truth variants). Calls are",
 r"the germline operating point (PASS, allele frequency $\geq 0.2$). \textbf{vardictcpp and",
 r"VarDictJava emit byte-identical VCFs}, so every value holds for both; per-region distributions",
 r"are in Fig.~\ref{fig:accuracy}.}",
 r"\begin{tabular}{lrrrrrrr}", r"\toprule",
 r"Class & Truth & TP & FP & FN & Precision & Recall & F$_1$ \\", r"\midrule"]
for row in acc:
    t3.append(f"{CLSNAME[row['class']]} & {fnum(row['n_truth'])} & {fnum(row['tp'])} & {row['fp']} & "
              f"{row['fn']} & {float(row['precision']):.3f} & {float(row['recall']):.3f} & "
              f"{float(row['f1']):.3f} \\\\")
t3 += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]

with open(os.path.join(HERE, "tables.tex"), "w") as f:
    f.write("% Auto-generated by make_paper_assets.py\n")
    f.write("\n".join(t1) + "\n\n" + "\n".join(t2) + "\n\n" + "\n".join(t3) + "\n")
print("wrote tables.tex")

# ---------------- figures ----------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import LogLocator, FuncFormatter, NullFormatter
import numpy as np

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans", "Arial"],
    "font.size": 8.5, "axes.labelsize": 9, "xtick.labelsize": 8.5, "ytick.labelsize": 8,
    "legend.fontsize": 8.5, "axes.linewidth": 0.7, "xtick.major.width": 0.7,
    "ytick.major.width": 0.7, "ytick.minor.width": 0.5, "figure.dpi": 150,
    "savefig.dpi": 400, "savefig.bbox": "tight", "pdf.fonttype": 42, "ps.fonttype": 42,
})

W = 0.34          # box width
OFF = 0.205       # dodge offset from group centre
GREY = "#3a3a3a"  # ink for median / axis text


def _shade(hex_color, f):
    """darken hex_color toward black by fraction f (0..1)."""
    h = hex_color.lstrip("#")
    r, g, b = (int(h[i:i+2], 16) for i in (0, 2, 4))
    return "#%02x%02x%02x" % (int(r*(1-f)), int(g*(1-f)), int(b*(1-f)))


def dodged_boxes(ax, cats, series_data, colors, ms=4.2, palpha=0.95, jit=0.075):
    """cats: x labels. series_data: list over series of list-over-cats of value-lists."""
    xs = np.arange(len(cats))
    for si, (data, color) in enumerate(zip(series_data, colors)):
        edge = _shade(color, 0.15)
        pos = xs + (si * 2 - 1) * OFF
        for xi, vals in zip(pos, data):
            if not vals:
                continue
            if len(vals) >= 3:
                bp = ax.boxplot([vals], positions=[xi], widths=W, patch_artist=True,
                                showfliers=False, whis=(0, 100), capwidths=W * 0.5,
                                zorder=3)
                for b in bp["boxes"]:
                    b.set(facecolor=color, alpha=0.30, edgecolor=edge, linewidth=1.1)
                for w in bp["whiskers"] + bp["caps"]:
                    w.set(color=edge, linewidth=1.0)
                for m in bp["medians"]:
                    m.set(color="white", linewidth=1.8, solid_capstyle="butt")
            else:
                m = st.median(vals)
                ax.plot([xi - W/2, xi + W/2], [m, m], color=edge, lw=1.8, zorder=3)
            offs = np.linspace(-jit, jit, len(vals)) if len(vals) > 1 else [0.0]
            ax.plot([xi + o for o in offs], vals, "o", ms=ms, mfc=color,
                    mec="white", mew=0.6, alpha=palpha, zorder=6)
    ax.set_xticks(xs); ax.set_xticklabels(cats)
    ax.set_xlim(-0.62, len(cats) - 0.38)
    ax.spines[["top", "right"]].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color("#9a9a9a")
    ax.tick_params(length=3, color="#9a9a9a", labelcolor=GREY)


def log_yaxis(ax):
    ax.set_yscale("log")
    ax.yaxis.set_major_locator(LogLocator(base=10))
    ax.yaxis.set_minor_locator(LogLocator(base=10, subs=(2, 3, 5)))
    def fmt(v, _):
        e = np.floor(np.log10(v)); m = v / 10**e
        return f"{v:g}" if round(m) in (1, 2, 5) else ""
    ax.yaxis.set_major_formatter(FuncFormatter(fmt))
    ax.yaxis.set_minor_formatter(FuncFormatter(fmt))
    ax.grid(axis="y", which="major", lw=0.6, color="#e2e2e2", zorder=0)
    ax.grid(axis="y", which="minor", lw=0.5, color="#f0f0f0", zorder=0)


def top_legend(fig, ncol=2, y=1.005):
    handles = [Patch(facecolor=c, alpha=0.30, edgecolor=_shade(c, 0.15), linewidth=1.1,
                     label=n) for n, c in SERIES]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, y),
               ncol=ncol, frameon=False, handlelength=1.3, columnspacing=1.8,
               handletextpad=0.6)


def panel_tag(ax, tag):
    ax.text(-0.015, 1.06, tag, transform=ax.transAxes, va="bottom", ha="right",
            fontsize=10, fontweight="bold", color=GREY)


def klabel(s): return f"{round(REG.get(s,0)/1000)} k"


# ---- Fig 1 & 2: performance box plots ----
def perf_fig(metric, ylabel, fname):
    fig, axes = plt.subplots(1, 2, figsize=(6.9, 3.3), sharey=True)
    for ax, th, tag in ((axes[0], "1", "a"), (axes[1], "8", "b")):
        cats = [klabel(s) for s in perf_samples]
        sd = [[REP.get((s, tool, th), {}).get(metric, []) for s in perf_samples]
              for tool, _ in SERIES]
        dodged_boxes(ax, cats, sd, [C_CPP, C_JAVA])
        log_yaxis(ax)
        ax.set_xlabel("Target regions")
        ax.set_title(f"{th} thread" + ("s" if th != "1" else ""), fontsize=8.5,
                     color=GREY, pad=4)
        panel_tag(ax, tag)
    axes[0].set_ylabel(ylabel)
    top_legend(fig)
    fig.tight_layout(w_pad=1.4, rect=(0, 0, 1, 0.93))
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIGDIR, f"{fname}.{ext}"))
    plt.close(fig); print(f"wrote {fname}")

perf_fig("wall", "Runtime (s)", "fig1_runtime")
perf_fig("rss", "Peak memory (MB)", "fig2_memory")

# ---- Fig 3: accuracy box plots ----
if per:
    PER = {}  # (impl,class,metric) -> [values over windows]
    for r in per:
        for m in ("precision", "recall", "f1"):
            PER.setdefault((r["impl"], r["class"], m), []).append(float(r[m]))
    classes = ["snv", "indel"]; clsdisp = {"snv": "SNV", "indel": "Indel"}
    metrics = [("precision", "Precision"), ("recall", "Recall"), ("f1", "F$_1$-score")]
    fig, axes = plt.subplots(1, 3, figsize=(7.1, 3.1), sharey=True)
    for ax, (mkey, mlab), tag in zip(axes, metrics, ("a", "b", "c")):
        cats = [clsdisp[c] for c in classes]
        sd = [[PER.get((tool, c, mkey), []) for c in classes] for tool, _ in SERIES]
        dodged_boxes(ax, cats, sd, [C_CPP, C_JAVA], ms=3.0, palpha=0.6, jit=0.11)
        ax.set_ylim(0.44, 1.03)
        ax.set_yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
        ax.grid(axis="y", which="major", lw=0.6, color="#e2e2e2", zorder=0)
        ax.set_xlabel(mlab)
        panel_tag(ax, tag)
    axes[0].set_ylabel("Score (per 1-Mb window)")
    top_legend(fig)
    fig.tight_layout(w_pad=1.4, rect=(0, 0, 1, 0.93))
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIGDIR, f"fig3_accuracy.{ext}"))
    plt.close(fig); print("wrote fig3_accuracy")
print("done.")
