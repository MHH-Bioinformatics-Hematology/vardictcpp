# attic — verified-but-unwired code

## cigar_modifier.{hpp,cpp}
Faithful port of VarDict's CigarModifier (read-end mismatch -> soft-clip, chimeric-seed clip removal).
VERIFIED CORRECT against instrumented VarDict: on the 1Mb/300x synthetic it reproduces 27601/27606 of
VarDict's CIGAR reshapings byte-for-byte (the 5 misses are indel reads, whose indel-collapse loop is
not ported).

NOT wired into the build because CigarModifier alone regresses counting (synthetic AltDepth 99.7%->90.7%):
its soft-clipping removes read-end bases that VarDict then MERGES BACK via realignment of the resulting
short (~3 bp) soft-clips. That merge-back is a coupled step not yet ported (the short clips fall below the
realignlgins/realignlgdel length thresholds). To enable CigarModifier, port the soft-clip-consensus
merge-back together with it, then re-verify counting stays >= 99.7%. This drives the ExtraAF/QStd columns.
