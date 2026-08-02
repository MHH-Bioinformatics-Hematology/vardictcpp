#!/usr/bin/env python3
"""Per-region precision/recall/F1 from rtg vcfeval output, for box plots.

Bins chr20 into fixed windows and, per window, computes:
  precision = TP_call / (TP_call + FP)      from tp.vcf.gz + fp.vcf.gz
  recall    = TP_base / (TP_base + FN)      from tp-baseline.vcf.gz + fn.vcf.gz
  F1        = 2PR/(P+R)
for each implementation (vardictcpp, VarDictJava) and variant class (all/SNV/indel).
Windows are kept only when both denominators have >= MIN_N variants, so single-variant
windows don't inject 0/1 noise. Emits a long-format TSV to stdout.
"""
import gzip, sys, os

WIN = 1_000_000
MIN_N = 5
HERE = os.path.dirname(os.path.abspath(__file__))

def vclass(ref, alt):
    return "snv" if (len(ref) == 1 and len(alt) == 1) else "indel"

def read(path):
    # yield (pos, class) for each record
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            c = ln.split("\t")
            pos = int(c[1]); ref = c[3]; alt = c[4].split(",")[0]
            yield pos, vclass(ref, alt)

def tally(path):
    # window -> {all,snv,indel} counts
    d = {}
    for pos, cls in read(path):
        w = pos // WIN
        b = d.setdefault(w, {"all": 0, "snv": 0, "indel": 0})
        b["all"] += 1; b[cls] += 1
    return d

def f1(p, r):
    return 0.0 if (p + r) == 0 else 2 * p * r / (p + r)

def main():
    impls = {"vardictcpp": "eval_cpp_germ", "VarDictJava": "eval_java_germ"}
    print("impl\tclass\twindow\tprecision\trecall\tf1\tn_call\tn_truth")
    for impl, d in impls.items():
        tp_c = tally(os.path.join(HERE, d, "tp.vcf.gz"))
        fp = tally(os.path.join(HERE, d, "fp.vcf.gz"))
        tp_b = tally(os.path.join(HERE, d, "tp-baseline.vcf.gz"))
        fn = tally(os.path.join(HERE, d, "fn.vcf.gz"))
        wins = set(tp_c) | set(fp) | set(tp_b) | set(fn)
        for w in sorted(wins):
            for cls in ("all", "snv", "indel"):
                g = lambda t: t.get(w, {}).get(cls, 0)
                tpc, fpp, tpb, fnn = g(tp_c), g(fp), g(tp_b), g(fn)
                ncall = tpc + fpp; ntruth = tpb + fnn
                if ncall < MIN_N or ntruth < MIN_N:
                    continue
                p = tpc / ncall
                r = tpb / ntruth
                print(f"{impl}\t{cls}\t{w}\t{p:.4f}\t{r:.4f}\t{f1(p, r):.4f}\t{ncall}\t{ntruth}")

if __name__ == "__main__":
    main()
