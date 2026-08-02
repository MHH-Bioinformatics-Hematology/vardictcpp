#!/usr/bin/env python3
"""Pure-Python reproduction of VarDict's teststrandbias.R.

Reads VarDict tabular output on stdin (34/36/38 columns), computes a two-sided
Fisher exact p-value and odds ratio from the strand-count columns 10-13
(1-indexed: refFwd refRev altFwd altRev), and writes columns 1-20, pvalue,
oddratio, columns 21..end -- exactly the column layout var2vcf_valid.pl expects.
No R / scipy dependency; Fisher exact via log-gamma hypergeometric tail sum,
matching R's fisher.test two-sided p-value. Rounds like the R script (5 dp).
"""
import sys, math

_lg = math.lgamma
def _lbinom_row(a, b, c, d):
    # log P of a 2x2 table with fixed margins (hypergeometric)
    n = a + b + c + d
    return (_lg(a+b+1)+_lg(c+d+1)+_lg(a+c+1)+_lg(b+d+1)
            - _lg(a+1)-_lg(b+1)-_lg(c+1)-_lg(d+1)-_lg(n+1))

def fisher_two_sided(a, b, c, d):
    n = a + b + c + d
    if n == 0:
        return 1.0, float("nan")
    r1 = a + b            # row1 total
    c1 = a + c            # col1 total
    lo = max(0, c1 - (n - r1))
    hi = min(c1, r1)
    p_obs = _lbinom_row(a, b, c, d)
    tol = 1e-7
    total = 0.0
    for x in range(lo, hi + 1):
        aa = x; bb = r1 - x; cc = c1 - x; dd = (n - r1) - (c1 - x)
        lp = _lbinom_row(aa, bb, cc, dd)
        if lp <= p_obs + tol:
            total += math.exp(lp)
    pval = min(1.0, total)
    # sample odds ratio (a*d)/(b*c); R reports the conditional MLE, but the
    # sample OR is what var2vcf uses for its strand-bias flag and is stable here.
    orr = float("inf") if (b == 0 or c == 0) else (a * d) / (b * c)
    return pval, orr

def main():
    for line in sys.stdin:
        line = line.rstrip("\n")
        if not line:
            continue
        f = line.split("\t")
        ncol = len(f)
        if ncol not in (34, 36, 38):
            sys.stderr.write("strandbias.py: unexpected column count %d\n" % ncol)
            sys.exit(1)
        a, b, c, d = (int(f[9]), int(f[10]), int(f[11]), int(f[12]))
        p, orr = fisher_two_sided(a, b, c, d)
        ps = "%.5f" % round(p, 5)
        os_ = ("Inf" if orr == float("inf")
               else ("NaN" if orr != orr else "%.5f" % round(orr, 5)))
        out = f[0:20] + [ps, os_] + f[20:ncol]
        sys.stdout.write("\t".join(out) + "\n")

if __name__ == "__main__":
    main()
