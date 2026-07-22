#!/bin/bash
# vardictcpp test suite: run against a small self-contained fixture and assert parity with the golden
# output produced by stock VarDictJava 1.8.3 (test/data/golden.simple.tsv).
#
# Asserts (all must hold):
#   1. Variant SET identical to VarDict  (0 false-positives, 0 false-negatives on chr/pos/ref/alt)
#   2. Depth and AltDepth identical on every shared variant
#   3. Full-row byte-identity >= 95% of VarDict's rows
set -u
DIR="$(cd "$(dirname "$0")" && pwd)"
VC="${VC:-$DIR/../build/vardictcpp}"
REF="$DIR/data/tref.fa"; BAM="$DIR/data/tref.bam"; GOLD="$DIR/data/golden.simple.tsv"
[ -x "$VC" ] || { echo "FAIL: vardictcpp binary not found at $VC (build first)"; exit 1; }

OUT="$(mktemp)"
"$VC" -G "$REF" -f 0.05 -N t -b "$BAM" -R tref:1-20001 -th 1 2>/dev/null | sort > "$OUT"

fail=0
# 1. variant set
setFP=$(comm -13 <(cut -f3,4,6,7 "$GOLD"|sort -u) <(cut -f3,4,6,7 "$OUT"|sort -u) | wc -l)
setFN=$(comm -23 <(cut -f3,4,6,7 "$GOLD"|sort -u) <(cut -f3,4,6,7 "$OUT"|sort -u) | wc -l)
[ "$setFP" -eq 0 ] && [ "$setFN" -eq 0 ] && echo "PASS variant-set: 0 FP, 0 FN" || { echo "FAIL variant-set: $setFP FP, $setFN FN"; fail=1; }

# 2. Depth (col8) + AltDepth (col9) on shared (chr,pos,ref,alt)
mism=$(python3 - "$GOLD" "$OUT" <<'PY'
import sys
g={}
for l in open(sys.argv[1]):
    f=l.rstrip('\n').split('\t'); g[(f[2],f[3],f[5],f[6])]=(f[7],f[8])
bad=0
for l in open(sys.argv[2]):
    f=l.rstrip('\n').split('\t'); k=(f[2],f[3],f[5],f[6])
    if k in g and (f[7],f[8])!=g[k]: bad+=1
print(bad)
PY
)
[ "$mism" -eq 0 ] && echo "PASS Depth/AltDepth: exact on all shared variants" || { echo "FAIL Depth/AltDepth: $mism mismatches"; fail=1; }

# 3. byte-identical rows
G=$(wc -l < "$GOLD"); BI=$(comm -12 "$GOLD" "$OUT" | wc -l)
pct=$(awk "BEGIN{printf \"%.1f\", 100*$BI/$G}")
awk "BEGIN{exit !($BI*100 >= $G*95)}" && echo "PASS byte-identity: $BI/$G ($pct%) >= 95%" || { echo "FAIL byte-identity: $pct% < 95%"; fail=1; }

rm -f "$OUT"
[ "$fail" -eq 0 ] && echo "== ALL TESTS PASSED ==" || echo "== TESTS FAILED =="
exit $fail
