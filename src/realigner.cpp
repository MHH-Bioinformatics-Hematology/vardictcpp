#include "realigner.hpp"
#include "util.hpp"
#include <algorithm>
#include <vector>
#include <string>
#include <set>
#include <cmath>

namespace vardict {

// ---- small helpers (VariationUtils) -------------------------------------------------------------

static void correctCnt(Variation& v) {
    if (v.varsCount < 0) v.varsCount = 0;
    if (v.highQualityReadsCount < 0) v.highQualityReadsCount = 0;
    if (v.lowQualityReadsCount < 0) v.lowQualityReadsCount = 0;
    if (v.meanPosition < 0) v.meanPosition = 0;
    if (v.meanQuality < 0) v.meanQuality = 0;
    if (v.meanMappingQuality < 0) v.meanMappingQuality = 0;
    if (v.getDir(true) < 0) v.addDir(true, -v.getDir(true));
    if (v.getDir(false) < 0) v.addDir(false, -v.getDir(false));
}

// VariationUtils.adjCnt (3-arg). referenceVar may be null.
static void adjCnt(Variation& varToAdd, const Variation& variant, Variation* referenceVar) {
    varToAdd.varsCount += variant.varsCount;
    varToAdd.extracnt += variant.varsCount;
    varToAdd.highQualityReadsCount += variant.highQualityReadsCount;
    varToAdd.lowQualityReadsCount += variant.lowQualityReadsCount;
    varToAdd.meanPosition += variant.meanPosition;
    varToAdd.meanQuality += variant.meanQuality;
    varToAdd.meanMappingQuality += variant.meanMappingQuality;
    varToAdd.numberOfMismatches += variant.numberOfMismatches;
    varToAdd.pstd = true;
    varToAdd.qstd = true;
    varToAdd.addDir(true, variant.getDir(true));
    varToAdd.addDir(false, variant.getDir(false));
    if (!referenceVar) return;
    referenceVar->varsCount -= variant.varsCount;
    referenceVar->highQualityReadsCount -= variant.highQualityReadsCount;
    referenceVar->lowQualityReadsCount -= variant.lowQualityReadsCount;
    referenceVar->meanPosition -= variant.meanPosition;
    referenceVar->meanQuality -= variant.meanQuality;
    referenceVar->meanMappingQuality -= variant.meanMappingQuality;
    referenceVar->numberOfMismatches -= variant.numberOfMismatches;
    referenceVar->subDir(true, variant.getDir(true));
    referenceVar->subDir(false, variant.getDir(false));
    correctCnt(*referenceVar);
}
static void adjCnt(Variation& a, const Variation& b) { adjCnt(a, b, nullptr); }

// VariationRealigner.adjRefCnt
static void adjRefCnt(Variation& tv, Variation* ref, int len) {
    if (!ref) return;
    double mp = tv.varsCount ? tv.meanPosition / tv.varsCount : 0;
    double f = tv.meanPosition != 0 ? (mp - len + 1) / mp : 0;
    if (f < 0) return;
    if (f > 1) f = 1;
    ref->varsCount -= (int)(f * tv.varsCount);
    ref->highQualityReadsCount -= (int)(f * tv.highQualityReadsCount);
    ref->lowQualityReadsCount -= (int)(f * tv.lowQualityReadsCount);
    ref->meanPosition -= f * tv.meanPosition;
    ref->meanQuality -= f * tv.meanQuality;
    ref->meanMappingQuality -= f * tv.meanMappingQuality;
    ref->numberOfMismatches -= f * tv.numberOfMismatches;
    ref->subDir(true, (int)(f * tv.getDir(true)));
    ref->subDir(false, (int)(f * tv.getDir(false)));
    correctCnt(*ref);
}

// VariationRealigner.adjRefFactor
static void adjRefFactor(Variation* ref, double f) {
    if (!ref) return;
    if (f > 1) f = 1;
    if (f < -1) return;
    int oldCnt = ref->varsCount;
    ref->varsCount -= (int)(f * ref->varsCount);
    ref->highQualityReadsCount -= (int)(f * ref->highQualityReadsCount);
    ref->lowQualityReadsCount -= (int)(f * ref->lowQualityReadsCount);
    double factorCnt = oldCnt != 0 ? std::abs(ref->varsCount - oldCnt) / (double)oldCnt : 1;
    if ((f < 0 && factorCnt > 0) || (f > 0 && factorCnt < 0)) factorCnt = -factorCnt;
    ref->meanPosition -= ref->meanPosition * factorCnt;
    ref->meanQuality -= ref->meanQuality * factorCnt;
    ref->meanMappingQuality -= ref->meanMappingQuality * factorCnt;
    ref->numberOfMismatches -= f * ref->numberOfMismatches;
    ref->varsCountOnForward -= (int)(f * ref->varsCountOnForward);
    ref->varsCountOnReverse -= (int)(f * ref->varsCountOnReverse);
    correctCnt(*ref);
}

static std::string joinRef(Reference& ref, int from, int to) {
    std::string s;
    for (int i = from; i <= to; ++i) if (ref.has(i)) s += ref.at(i);
    return s;
}

// VariationUtils.joinRef(Map,int,double) overload: iterates with a STRICT `i < to` (exclusive upper
// bound), unlike the int overload's inclusive `i <= to`. Used by realignlgins30's tandem-duplication
// branch, where the bound is a floating-point expression.
static std::string joinRef(Reference& ref, int from, double to) {
    std::string s;
    for (int i = from; i < to; ++i) if (ref.has(i)) s += ref.at(i);
    return s;
}
static std::string joinRefFor5Lgins(Reference& ref, int from, int to, const std::string& seq, const std::string& EXTRA) {
    std::string sb;
    for (int i = from; i <= to; ++i) {
        if (to - i < (int)seq.size() - (int)EXTRA.size()) { char c = charAt(seq, to - i + (int)EXTRA.size()); if (c != (char)-1) sb += c; }
        else if (ref.has(i)) sb += ref.at(i);
    }
    return sb;
}
static std::string joinRefFor3Lgins(Reference& ref, int from, int to, int shift5, const std::string& seq, const std::string& EXTRA) {
    std::string sb;
    for (int i = from; i <= to; ++i) {
        if (i - from >= shift5 && i - from - shift5 < (int)seq.size() - (int)EXTRA.size()) {
            char c = charAt(seq, i - from - shift5 + (int)EXTRA.size()); if (c != (char)-1) sb += c;
        } else if (ref.has(i)) sb += ref.at(i);
    }
    return sb;
}

static Variation& getVariation(std::unordered_map<int, VarMap>& hash, int pos, const std::string& key) {
    return hash[pos][key];
}
static Variation* getVariationMaybe(std::unordered_map<int, VarMap>& hash, int pos, char refBase) {
    auto it = hash.find(pos);
    if (it == hash.end()) return nullptr;
    auto vit = it->second.find(std::string(1, refBase));
    return vit == it->second.end() ? nullptr : &vit->second;
}

static int countChar(const std::string& s, char c) { return (int)std::count(s.begin(), s.end(), c); }
static bool islowcomplexseq(const std::string& seq) {
    int len = (int)seq.size();
    if (len == 0) return true;
    int ntcnt = 0;
    for (char c : {'A','T','G','C'}) {
        int k = countChar(seq, c);
        if (k > 0) ntcnt++;
        if (k / (double)len > 0.75) return true;
    }
    return ntcnt < 3;
}

// ismatch: two sequences match with <=MM mismatches and <15% mismatch fraction. dir is +1 or -1.
static bool ismatch(std::string seq1, std::string seq2, int dir, int MM = 3) {
    std::string s2;
    for (char c : seq2) if (c != '#' && c != '^') s2 += c;   // strip #,^
    int mm = 0;
    for (int n = 0; n < (int)seq1.size() && n < (int)s2.size(); ++n) {
        char c2 = charAt(s2, dir * n - (dir == -1 ? 1 : 0));
        if (seq1[n] != c2) mm++;
    }
    return mm <= MM && mm / (double)seq1.size() < 0.15;
}

// VariationRealigner.rmCnt: subtract tv's tallies from vref.
static void rmCnt(Variation& vref, const Variation& tv) {
    vref.varsCount -= tv.varsCount;
    vref.highQualityReadsCount -= tv.highQualityReadsCount;
    vref.lowQualityReadsCount -= tv.lowQualityReadsCount;
    vref.meanPosition -= tv.meanPosition;
    vref.meanQuality -= tv.meanQuality;
    vref.meanMappingQuality -= tv.meanMappingQuality;
    vref.subDir(true, tv.getDir(true));
    vref.subDir(false, tv.getDir(false));
    correctCnt(vref);
}

// VariationRealigner.findbp: slide `sequence` against the reference within indelsize to find a
// breakpoint where it matches with <=3 mismatches (large-indel breakpoint detection).
static int findbp(const std::string& sequence, int startPosition, Reference& ref, int direction,
                  int chrLen, int indelsize) {
    const int maxmm = 3;
    int bp = 0, score = 0;
    int idx = chrLen;
    for (int n = 0; n < indelsize; ++n) {
        int mm = 0, i = 0;
        std::set<char> m;
        for (i = 0; i < (int)sequence.size(); ++i) {
            int rp = startPosition + direction * n + direction * i;
            if (rp < 1 || rp > idx) break;
            if (ref.has(rp) && sequence[i] == ref.at(rp)) m.insert(sequence[i]);
            else mm++;
            if (mm > maxmm - n / 100) break;
        }
        if ((int)m.size() < 3) continue;
        if (mm <= maxmm - n / 100 && i >= (int)sequence.size() - 2 && i >= 8 + n / 10 &&
            mm / (double)i < 0.12) {
            int lbp = startPosition + direction * n - (direction < 0 ? direction : 0);
            if (mm == 0 && i == (int)sequence.size()) return lbp;
            else if (i - mm > score) { bp = lbp; score = i - mm; }
        }
    }
    return bp;
}

// VariationRealigner.ismatchref: does `sequence` match the reference at `position` (stepping by dir)
// with <=MM mismatches and <15% mismatch fraction?
static bool ismatchref(const std::string& sequence, Reference& ref, int position, int dir, int MM = 3) {
    int mm = 0;
    for (int n = 0; n < (int)sequence.size(); ++n) {
        int rp = position + dir * n;
        if (!ref.has(rp)) return false;
        char sc = charAt(sequence, dir == 1 ? n : dir * n - 1);
        if (sc != ref.at(rp)) mm++;
    }
    return mm <= MM && mm / (double)sequence.size() < 0.15;
}

// StructuralVariantsProcessor.findMatch: seed-lookup the soft-clip consensus in the reference to
// locate a breakpoint elsewhere in the loaded window (used when findbp fails). Returns {bp, extra}.
struct Match { int bp; std::string extra; };
static Match findMatch(std::string seq, Reference& ref, int /*position*/, int dir, int SEED, int MM) {
    if (dir == -1) seq = reverseStr(seq);
    std::string extra;
    for (int i = (int)seq.size() - SEED; i >= 0; --i) {
        std::string kmer = substr(seq, i, SEED);
        int firstSeed = ref.seedUnique(kmer);
        if (firstSeed <= 0) continue;
        int bp = dir == 1 ? firstSeed - i : firstSeed + (int)seq.size() - i - 1;
        if (ismatchref(seq, ref, bp, dir, MM)) {
            int mm = dir == -1 ? -1 : 0;
            while (ref.has(bp) && charAt(seq, mm) != (char)-1 && ref.at(bp) != charAt(seq, mm)) {
                extra += substr(seq, mm, 1); bp += dir; mm += dir;
            }
            if (!extra.empty() && dir == -1) extra = reverseStr(extra);
            return { bp, extra };
        } else {
            // Complex-indel fallback: walk up to 15bp from the seed, allowing mismatched end bases
            // (collected into EXTRA), until the remaining consensus matches the reference (MM=1).
            auto hasNe = [&](char c, int pos) { return ref.has(pos) && ref.at(pos) != c; };
            std::string sseq = seq;
            int eqcnt = 0;
            for (int ii = 1; ii <= 15; ++ii) {
                bp += dir;
                sseq = dir == 1 ? substr(sseq, 1) : substr(sseq, 0, (int)sseq.size() - 1);
                if (dir == 1) {
                    if (hasNe(charAt(sseq, 0), bp)) continue;
                    eqcnt++;
                    if (hasNe(charAt(sseq, 1), bp + 1)) continue;
                    extra = substr(seq, 0, ii);
                } else {
                    if (hasNe(charAt(sseq, -1), bp)) continue;
                    eqcnt++;
                    if (hasNe(charAt(sseq, -2), bp - 1)) continue;
                    extra = substr(seq, -ii);
                }
                if (eqcnt >= 3 && eqcnt / (double)ii > 0.5) break;
                if (ismatchref(sseq, ref, bp, dir, 1)) return { bp, extra };
            }
        }
    }
    return { 0, "" };
}

// StructuralVariantsProcessor.findMatchRev: like findMatch but searches the reverse-complemented
// consensus (used by findsv's candidate-inversion path). dir==1 means the seq came from a 3' clip.
static Match findMatchRev(std::string seq, Reference& ref, int /*position*/, int dir, int SEED, int MM) {
    if (dir == 1) seq = reverseStr(seq);
    seq = complementStr(seq);
    auto hasNe = [&](char c, int pos) { return ref.has(pos) && ref.at(pos) != c; };
    std::string extra;
    for (int i = (int)seq.size() - SEED; i >= 0; --i) {
        std::string kmer = substr(seq, i, SEED);
        int firstSeed = ref.seedUnique(kmer);
        if (firstSeed <= 0) continue;
        int bp = dir == 1 ? firstSeed + (int)seq.size() - i - 1 : firstSeed - i;
        if (ismatchref(seq, ref, bp, -dir, MM)) {
            return { bp, extra };
        } else {
            // for complex indels, allow some mismatches at the end up to 15bp or 20% length
            std::string sseq = seq;
            int eqcnt = 0;
            for (int j = 1; j <= 15; ++j) {
                bp -= dir;
                sseq = dir == -1 ? substr(sseq, 1) : substr(sseq, 0, (int)sseq.size() - 1);
                if (dir == -1) {
                    if (hasNe(charAt(sseq, 0), bp)) continue;
                    eqcnt++;
                    if (hasNe(charAt(sseq, 1), bp + 1)) continue;
                    extra = substr(seq, 0, j);
                } else {
                    if (hasNe(charAt(sseq, -1), bp)) continue;
                    eqcnt++;
                    if (hasNe(charAt(sseq, -2), bp - 1)) continue;
                    extra = substr(seq, -j);
                }
                if (eqcnt >= 3 && eqcnt / (double)j > 0.5) break;
                if (ismatchref(sseq, ref, bp, -dir, 1)) return { bp, extra };
            }
        }
    }
    return { 0, "" };
}
static Match findMatchRev(std::string seq, Reference& ref, int position, int dir) {
    return findMatchRev(std::move(seq), ref, position, dir, Reference::SEED_1, 3);
}

// ---- soft-clip consensus (findconseq) -----------------------------------------------------------

static std::string findconseq(Sclip& sc) {
    if (!sc.sequence.empty() || sc.used) return sc.sequence; // sequence cached (may be "")
    int total = 0, match = 0;
    std::string seq;
    bool flag = false;
    for (auto& [posInSclip, nv] : sc.nt) {
        int maxCount = 0; double maxQuality = 0; char chosen = 0; int totalCount = 0;
        for (auto& [base, cnt] : nv) {
            totalCount += cnt;
            double mq = 0;
            auto sit = sc.seq.find(posInSclip);
            bool hasQ = sit != sc.seq.end() && sit->second.count(base);
            if (hasQ) mq = sit->second[base]->meanQuality;
            if (cnt > maxCount || (hasQ && mq > maxQuality)) {
                maxCount = cnt; chosen = base; maxQuality = mq;
            }
        }
        if (posInSclip == 3 && (int)sc.nt.size() >= 6 &&
            totalCount / (double)sc.varsCount < 0.2 && totalCount <= 2) break;
        if ((totalCount - maxCount > 2 || maxCount <= totalCount - maxCount) &&
            maxCount / (double)totalCount < 0.8) {
            if (flag) break;
            flag = true;
        }
        total += totalCount; match += maxCount;
        if (chosen) seq += chosen;
    }
    int ntSize = (int)sc.nt.size();
    std::string SEQ;
    if (total != 0 && match / (double)total > 0.9 &&
        seq.size() / 1.5 > ntSize - (int)seq.size() &&
        (seq.size() / (double)ntSize > 0.8 || ntSize - (int)seq.size() < 10 || seq.size() > 25)) {
        SEQ = seq;
    }
    if (!SEQ.empty() && (int)SEQ.size() > Reference::SEED_2) {
        // B_A7 = ^.AAAAAAA, B_T7 = ^.TTTTTTT: poly-A/T runs near the clip start -> unusable consensus.
        if (SEQ.size() >= 8 && (SEQ.substr(1, 7) == "AAAAAAA" || SEQ.substr(1, 7) == "TTTTTTT")) sc.used = true;
        if (islowcomplexseq(SEQ)) sc.used = true;
    }
    sc.sequence = SEQ;
    return SEQ;
}

// ---- mismatch scanning (findMM5 / findMM3) -------------------------------------------------------

struct Mismatch { std::string seq; int pos; int end; };
struct MismatchResult { std::vector<Mismatch> mm; std::vector<int> scp; int nm; int misp; std::string misnt; };

static MismatchResult findMM5(VariationData& vd, Reference& ref, int position, std::string wupseq) {
    std::string seq;
    for (char c : wupseq) if (c != '#' && c != '^') seq += c;
    const int longmm = 3;
    MismatchResult r; r.nm = 0; r.misp = 0;
    int n = 0, mn = 0, mcnt = 0;
    std::string str;
    while (ref.has(position - n) && charAt(seq, -1 - n) != (char)-1 &&
           ref.at(position - n) != charAt(seq, -1 - n) && mcnt < longmm) {
        str.insert(str.begin(), charAt(seq, -1 - n));
        r.mm.push_back({str, position - n, 5});
        n++; mcnt++;
    }
    r.scp.push_back(position + 1);
    char misnt = 0;
    if (str.size() == 1) {
        while (ref.has(position - n) && ref.at(position - n) == charAt(seq, -1 - n)) {
            n++; if (n != 0) mn++;
        }
        if (mn > 1) {
            int n2 = 0;
            while (-1 - n - 1 - n2 >= 0 && ref.has(position - n - 1 - n2) &&
                   ref.at(position - n - 1 - n2) == charAt(seq, -1 - n - 1 - n2)) n2++;
            if (n2 > 2) {
                r.scp.push_back(position - n - n2);
                r.misp = position - n; misnt = charAt(seq, -1 - n);
                auto it = vd.softClips5End.find(position - n - n2);
                if (it != vd.softClips5End.end()) it->second.used = true;
                mn += n2;
            } else {
                r.scp.push_back(position - n);
                auto it = vd.softClips5End.find(position - n);
                if (it != vd.softClips5End.end()) it->second.used = true;
            }
        }
    }
    r.nm = mn; if (misnt) r.misnt = std::string(1, misnt);
    return r;
}

static MismatchResult findMM3(VariationData& vd, Reference& ref, int p, std::string sanpseq) {
    std::string seq;
    for (char c : sanpseq) if (c != '#' && c != '^') seq += c;
    const int longmm = 3;
    MismatchResult r; r.nm = 0; r.misp = 0;
    int n = 0, mn = 0, mcnt = 0;
    std::string str;
    while (n < (int)seq.size() && ref.has(p + n) && ref.at(p + n) == seq[n]) n++;
    r.scp.push_back(p + n);
    int Tbp = p + n;
    while (mcnt <= longmm && n < (int)seq.size() && ref.has(p + n) && ref.at(p + n) != seq[n]) {
        str += seq[n];
        r.mm.push_back({str, Tbp, 3});
        n++; mcnt++;
    }
    char misnt = 0;
    if (str.size() == 1) {
        while (n < (int)seq.size() && ref.has(p + n) && ref.at(p + n) == seq[n]) { n++; if (n != 0) mn++; }
        if (mn > 1) {
            int n2 = 0;
            while (n + n2 + 1 < (int)seq.size() && ref.has(p + n + 1 + n2) && ref.at(p + n + 1 + n2) == seq[n + n2 + 1]) n2++;
            if (n2 > 2 && n + n2 + 1 < (int)seq.size()) {
                r.scp.push_back(p + n + n2);
                r.misp = p + n; misnt = seq[n];
                auto it = vd.softClips3End.find(p + n + n2);
                if (it != vd.softClips3End.end()) it->second.used = true;
                mn += n2;
            } else {
                r.scp.push_back(p + n);
                auto it = vd.softClips3End.find(p + n);
                if (it != vd.softClips3End.end()) it->second.used = true;
            }
        }
    }
    r.nm = mn; if (misnt) r.misnt = std::string(1, misnt);
    return r;
}

// ---- adjustMNP (unchanged) ----------------------------------------------------------------------

static std::string decorate(std::string s) { if (s.size() > 1) s.insert(1, "&"); return s; }

void adjustMNP(VariationData& vd, Reference& ref, const Config&, const Region&) {
    struct Item { int position; std::string desc; int count; };
    std::vector<Item> tmp;
    for (auto& [pos, m] : vd.mnp) for (auto& [desc, cnt] : m) tmp.push_back({pos, desc, cnt});
    std::sort(tmp.begin(), tmp.end(), [](const Item& a, const Item& b) {
        if (a.count != b.count) return a.count > b.count;
        if (a.position != b.position) return a.position < b.position;
        return a.desc > b.desc;
    });
    for (const auto& it : tmp) {
        int position = it.position; const std::string& vn = it.desc;
        auto pit = vd.nonInsertionVariants.find(position);
        if (pit == vd.nonInsertionVariants.end()) continue;
        VarMap& vp = pit->second;
        auto vrefIt = vp.find(vn);
        if (vrefIt == vp.end()) continue;
        Variation& vref = vrefIt->second;
        std::string mnt = vn; auto amp = mnt.find('&'); if (amp != std::string::npos) mnt.erase(amp, 1);
        for (int i = 0; i < (int)mnt.size() - 1; ++i) {
            std::string left = decorate(substr(mnt, 0, i + 1));
            std::string right = decorate(substr(mnt, -((int)mnt.size() - i - 1)));
            auto lit = vp.find(left);
            if (lit != vp.end()) {
                Variation& tref = lit->second;
                if (tref.varsCount > 0 && tref.varsCount < vref.varsCount &&
                    tref.meanPosition / tref.varsCount <= i + 1) { adjCnt(vref, tref); vp.erase(left); }
            }
            auto rposIt = vd.nonInsertionVariants.find(position + i + 1);
            if (rposIt != vd.nonInsertionVariants.end()) {
                auto rit = rposIt->second.find(right);
                if (rit != rposIt->second.end()) {
                    Variation& tref = rit->second;
                    if (tref.varsCount >= 0 && tref.varsCount < vref.varsCount) {
                        adjCnt(vref, tref); vd.refCoverage[position] += tref.varsCount; rposIt->second.erase(right);
                    }
                }
            }
        }
        // Absorb a 3' soft-clip consensus that starts with the MNP (the aligner clipped reads at the
        // MNP mismatches instead of aligning them). VariationRealigner.adjustMNP lines 346-358.
        {
            auto scIt = vd.softClips3End.find(position);
            if (scIt != vd.softClips3End.end() && !scIt->second.used) {
                Sclip& sc3v = scIt->second;
                std::string seq = findconseq(sc3v);
                if (seq.size() >= mnt.size() && seq.compare(0, mnt.size(), mnt) == 0) {
                    if (seq.size() == mnt.size() ||
                        ismatchref(seq.substr(mnt.size()), ref, position + (int)mnt.size(), 1)) {
                        adjCnt(vref, sc3v);
                        vd.refCoverage[position] += sc3v.varsCount;
                        sc3v.used = true;
                    }
                }
            }
        }
        // Absorb a 5' soft-clip consensus ending with the MNP (reverse orientation). Lines 360-377.
        {
            auto scIt = vd.softClips5End.find(position + (int)mnt.size());
            if (scIt != vd.softClips5End.end() && !scIt->second.used) {
                Sclip& sc5v = scIt->second;
                std::string seq = findconseq(sc5v);
                if (!seq.empty() && seq.size() >= mnt.size()) {
                    std::reverse(seq.begin(), seq.end());
                    if (seq.compare(seq.size() - mnt.size(), mnt.size(), mnt) == 0) {
                        if (seq.size() == mnt.size() ||
                            ismatchref(seq.substr(0, seq.size() - mnt.size()), ref, position - 1, -1)) {
                            adjCnt(vref, sc5v);
                            vd.refCoverage[position] += sc5v.varsCount;
                            sc5v.used = true;
                        }
                    }
                }
            }
        }
    }
}

// ---- realignins ---------------------------------------------------------------------------------

// One insertion allele of VariationRealigner.realignins: attracts nearby mismatched SNVs and matching
// soft clips into the "+SEQ" insertion. Extracted so realignlgins30 can re-run it on the single allele
// it just created (Java calls realignins({bi:{ins:cnt}}) inline), which the whole-map realignins pass
// -- run earlier in the pipeline -- never sees.
static void realignOneIns(VariationData& vd, Reference& ref, const Config& cfg, int maxReadLength,
                          int position, const std::string& vn, int insertionCount) {
        auto& NIV = vd.nonInsertionVariants;
        auto isATGC = [](char c){ return c=='A'||c=='C'||c=='G'||c=='T'; };
        // Parse the insertion-description grammar (VariationRealigner.realignins): a complex insertion
        // "+SEQ", optionally with "&extra", "#compm", "<dupN>ins3" or a "^N"/"^ATGC" tail. The old port
        // only handled plain "+SEQ" and dropped any '&'/'#' insertion, so complex insertions never
        // consumed their soft clips (and later plain insertions wrongly did).
        if (vn.empty() || vn[0] != '+') return;
        std::string insert;   // BEGIN_PLUS_ATGC: leading ATGC run after '+'
        for (size_t k = 1; k < vn.size() && isATGC(vn[k]); ++k) insert += vn[k];
        if (insert.empty()) return;
        int inslen = (int)insert.size();
        std::string ins3;     // DUP_NUM_ATGC: <dup(\d+)>([ATGC]+)$
        { size_t d = vn.find("<dup");
          if (d != std::string::npos) { size_t gt = vn.find('>', d);
            if (gt != std::string::npos) { int num = atoi(vn.c_str() + d + 4);
              std::string tail = vn.substr(gt + 1);
              bool allatgc = !tail.empty(); for (char c : tail) if (!isATGC(c)) allatgc = false;
              if (allatgc) { ins3 = tail; inslen += num + (int)ins3.size(); } } } }
        std::string extra;    // AMP_ATGC: &([ATGC]+)
        { size_t a = vn.find('&');
          if (a != std::string::npos) for (size_t k = a + 1; k < vn.size() && isATGC(vn[k]); ++k) extra += vn[k]; }
        std::string compm;    // HASH_ATGC: #([ATGC]+)
        { size_t h = vn.find('#');
          if (h != std::string::npos) for (size_t k = h + 1; k < vn.size() && isATGC(vn[k]); ++k) compm += vn[k]; }
        int newdel = 0;       // UP_NUMBER_END: \^(\d+)$
        { size_t c = vn.rfind('^');
          if (c != std::string::npos) { std::string tl = vn.substr(c + 1);
            bool alldig = !tl.empty(); for (char ch : tl) if (!isdigit((unsigned char)ch)) alldig = false;
            if (alldig) newdel = atoi(tl.c_str()); } }
        // tn = vn stripped of: leading '+', first '&', first '#', trailing "^\d+", first '^'.
        std::string tn = vn.substr(1);
        { size_t a = tn.find('&'); if (a != std::string::npos) tn.erase(a, 1); }
        { size_t h = tn.find('#'); if (h != std::string::npos) tn.erase(h, 1); }
        { size_t c = tn.rfind('^');
          if (c != std::string::npos) { bool alldig = c + 1 < tn.size();
            for (size_t k = c + 1; k < tn.size(); ++k) if (!isdigit((unsigned char)tn[k])) alldig = false;
            if (alldig) tn.erase(c); } }
        { size_t c = tn.find('^'); if (c != std::string::npos) tn.erase(c, 1); }

        int wustart = position - 150 > 1 ? position - 150 : 1;
        std::string wupseq = joinRef(ref, wustart, position) + tn;
        int sanend = position + (int)vn.size() + 100;
        std::string sanpseq;
        int mm3anchor;
        if (!ins3.empty()) {
            int p3 = position + inslen - (int)ins3.size() + Config::SVFLANK;
            if ((int)ins3.size() > Config::SVFLANK) sanpseq = substr(ins3, Config::SVFLANK - (int)ins3.size());
            sanpseq += joinRef(ref, position + 1, position + 101);
            mm3anchor = p3 + 1;
        } else {
            sanpseq = tn + joinRef(ref, position + (int)extra.size() + 1 + (int)compm.size() + newdel, sanend);
            mm3anchor = position + 1;
        }
        MismatchResult f3 = findMM3(vd, ref, mm3anchor, sanpseq);
        MismatchResult f5 = findMM5(vd, ref, position + (int)extra.size() + (int)compm.size() + newdel, wupseq);

        std::vector<Mismatch> mmm = f3.mm; mmm.insert(mmm.end(), f5.mm.begin(), f5.mm.end());
        Variation& vref = vd.insertionVariants[position][vn];
        for (auto& mismatch : mmm) {
            std::string mb = mismatch.seq; int mp = mismatch.pos; int me = mismatch.end;
            if (mb.size() > 1) mb = std::string(1, mb[0]) + "&" + mb.substr(1);
            auto pit = NIV.find(mp); if (pit == NIV.end()) continue;
            auto vit = pit->second.find(mb); if (vit == pit->second.end()) continue;
            Variation& variation = vit->second;
            if (variation.varsCount == 0) continue;
            if (variation.meanQuality / variation.varsCount < cfg.goodq) continue;
            if (variation.meanPosition / variation.varsCount > (me == 3 ? f3.nm + 4 : f5.nm + 4)) continue;
            if (variation.varsCount >= insertionCount + (int)insert.size() || variation.varsCount / insertionCount >= 8) continue;
            if (mp > position && me == 5) vd.refCoverage[position] += variation.varsCount;
            Variation* lref = nullptr;
            if (mp > position && me == 3 && ref.has(position))
                lref = getVariationMaybe(NIV, position, ref.at(position));
            adjCnt(vref, variation, lref);
            pit->second.erase(mb);
            if (pit->second.empty()) NIV.erase(pit);
        }
        if (f3.misp != 0 && f3.mm.size() == 1) {
            auto pit = NIV.find(f3.misp);
            if (pit != NIV.end()) { auto vit = pit->second.find(f3.misnt);
                if (vit != pit->second.end() && vit->second.varsCount < insertionCount) pit->second.erase(vit); }
        }
        if (f5.misp != 0 && f5.mm.size() == 1) {
            auto pit = NIV.find(f5.misp);
            if (pit != NIV.end()) { auto vit = pit->second.find(f5.misnt);
                if (vit != pit->second.end() && vit->second.varsCount < insertionCount) pit->second.erase(vit); }
        }
        for (int sc5pp : f5.scp) {
            auto it = vd.softClips5End.find(sc5pp);
            if (it == vd.softClips5End.end() || it->second.used) continue;
            Sclip& tv = it->second;
            std::string seq = findconseq(tv);
            if (!seq.empty() && ismatch(seq, wupseq, -1)) {
                if (sc5pp > position) vd.refCoverage[position] += tv.varsCount;
                adjCnt(vref, tv); tv.used = true;
            }
        }
        for (int sc3pp : f3.scp) {
            auto it = vd.softClips3End.find(sc3pp);
            if (it == vd.softClips3End.end() || it->second.used) continue;
            Sclip& tv = it->second;
            std::string seq = findconseq(tv);
            std::string mseq = !ins3.empty() ? sanpseq : substr(sanpseq, sc3pp - position - 1);
            if (!seq.empty() && ismatch(seq, mseq, 1)) {
                double tvmp = tv.varsCount ? tv.meanPosition / tv.varsCount : 0;
                if (sc3pp <= position || (int)insert.size() > tvmp) vd.refCoverage[position] += tv.varsCount;
                Variation* lref = nullptr;
                if (sc3pp > position && ref.has(position)) lref = getVariationMaybe(NIV, position, ref.at(position));
                if ((int)insert.size() > tvmp) lref = nullptr;
                adjCnt(vref, tv, lref); tv.used = true;
            }
        }
        if (!f3.scp.empty() && !f5.scp.empty()) {
            int first3 = f3.scp[0], first5 = f5.scp[0];
            if (first3 > first5 + 3 && first3 - first5 < maxReadLength * 0.75) {
                if (ref.has(position)) {
                    Variation* rv = getVariationMaybe(NIV, position, ref.at(position));
                    if (rv) adjRefFactor(rv, (first3 - first5 - 1) / (double)maxReadLength);
                }
                adjRefFactor(&vref, -(first3 - first5 - 1) / (double)maxReadLength);
            }
        }
}

void realignins(VariationData& vd, Reference& ref, const Config& cfg, const Region& region, int maxReadLength) {
    (void)region;
    struct Item { int position; std::string desc; int count; };
    std::vector<Item> tmp;
    for (auto& [pos, m] : vd.positionToInsertionCount) for (auto& [d, c] : m) tmp.push_back({pos, d, c});
    std::sort(tmp.begin(), tmp.end(), [](const Item& a, const Item& b) {
        if (a.count != b.count) return a.count > b.count;
        if (a.position != b.position) return a.position < b.position;
        return a.desc > b.desc;
    });
    for (const auto& t : tmp) realignOneIns(vd, ref, cfg, maxReadLength, t.position, t.desc, t.count);
    // Merge "+SEQ&extra" duplicates into "+SEQ" (ATGSs_AMP_ATGSs_END). Plain inserts skip this.
    for (int i = (int)tmp.size() - 1; i > 0; --i) {
        int p = tmp[i].position; const std::string& vn = tmp[i].desc;
        auto pit = vd.insertionVariants.find(p); if (pit == vd.insertionVariants.end()) continue;
        auto vit = pit->second.find(vn); if (vit == pit->second.end()) continue;
        auto amp = vn.find('&');
        if (amp != std::string::npos && vn[0] == '+') {
            std::string tnkey = vn.substr(0, amp);
            auto tit = pit->second.find(tnkey);
            if (tit != pit->second.end() && vit->second.varsCount < tit->second.varsCount) {
                Variation* lref = ref.has(p) ? getVariationMaybe(vd.nonInsertionVariants, p, ref.at(p)) : nullptr;
                adjCnt(tit->second, vit->second, lref); pit->second.erase(vit);
            }
        }
    }
}

// ---- realigndel ---------------------------------------------------------------------------------

// VariationRealigner.noPassingReads / vardict.pl sub noPassingReads. Scans the raw BAM over [s,e]
// (1-based inclusive, NO mapq/flag filtering, exactly like `samtools view $bam $chr:$s-$e`) and
// returns true when no read cleanly spans the gap as reference (cnt<=0). Reads whose CIGAR string
// already contains the exact "<e-s>D" deletion are skipped (they support the deletion). rlen is the
// aligned length (M+D only, CigarParser.getAlignedLength); re = alignmentStart + rlen. A "passing"
// read has re > e+2 && rs < s-2. Java's `midcnt` term (midcnt+1>0) is always true, so it is omitted.
static bool noPassingReads(const std::vector<BamReader*>& bams, const std::string& chr, int s, int e) {
    int cnt = 0;
    std::string dlenqr = std::to_string(e - s) + "D";
    for (BamReader* br : bams) {
        if (!br) continue;
        bam_hdr_t* hdr = br->hdr();
        int tid = bam_name2id(hdr, chr.c_str());
        if (tid < 0) {
            std::string alt = (chr.rfind("chr", 0) == 0) ? chr.substr(3) : "chr" + chr;
            tid = bam_name2id(hdr, alt.c_str());
        }
        if (tid < 0) continue;
        hts_itr_t* it = sam_itr_queryi(br->idx(), tid, s - 1, e);
        if (!it) continue;
        bam1_t* b = bam_init1();
        while (sam_itr_next(br->fp(), it, b) >= 0) {
            const bam1_core_t& c = b->core;
            if (c.n_cigar == 0) continue;
            const uint32_t* cig = bam_get_cigar(b);
            std::string cigstr;
            int rlen = 0;
            for (uint32_t k = 0; k < c.n_cigar; ++k) {
                int op = bam_cigar_op(cig[k]);
                int ol = bam_cigar_oplen(cig[k]);
                cigstr += std::to_string(ol);
                cigstr += bam_cigar_opchr(cig[k]);
                if (op == BAM_CMATCH || op == BAM_CDEL) rlen += ol;
            }
            if (cigstr.find(dlenqr) != std::string::npos) continue;
            int rs = (int)c.pos + 1; // 1-based alignment start
            int re = rs + rlen;
            if (re > e + 2 && rs < s - 2) cnt++;
        }
        bam_destroy1(b);
        hts_itr_destroy(it);
    }
    return cnt <= 0;
}

// One deletion's realignment: attribute nearby mismatch SNVs + soft-clip consensus to the deletion
// (VariationRealigner.realigndel per-position body). Extracted so realignlgdel can re-run it on the
// single large deletion it just built (Java realigndel(bams, dels5), l.1174/1345).
static void realignOneDel(VariationData& vd, Reference& ref, const Config& cfg, const Region& region,
                          int maxReadLength, const std::vector<BamReader*>& bams,
                          int p, const std::string& vn, int dcnt) {
    auto& NIV = vd.nonInsertionVariants;
    {
        // BEGIN_MINUS_NUMBER: dellen. Complex deletions "-N&ss"/"-N^ins"/"-N#seg^M" carry a tail that
        // shifts the 5'/3' flanking sequences used for soft-clip re-matching, so a complex deletion
        // does NOT scoop the soft-clip that belongs to the plain "-N" (VariationRealigner realigndel).
        //   BEGIN_MINUS_NUMBER_ANY  extra    = tail after "-<dellen>", with ^,&,# stripped
        //   CARET_ATGNC             extrains = ATGNC run following the first '^'
        //   UP_NUMBER_END           dellen  += trailing "^<digits>"
        if (vn.empty() || vn[0] != '-') return;
        int dellen = std::atoi(vn.c_str() + 1);
        std::string extra, extrains;
        {
            size_t di = 1; while (di < vn.size() && isdigit((unsigned char)vn[di])) di++;
            std::string rest = vn.substr(di);
            for (char ch : rest) if (ch != '^' && ch != '&' && ch != '#') extra += ch;
            auto cp = rest.find('^');
            if (cp != std::string::npos)
                for (size_t z = cp + 1; z < rest.size(); ++z) {
                    char c = rest[z];
                    if (c=='A'||c=='T'||c=='G'||c=='N'||c=='C') extrains += c; else break;
                }
            auto up = rest.rfind('^');
            if (up != std::string::npos && up + 1 < rest.size()) {
                bool allDig = true;
                for (size_t z = up + 1; z < rest.size(); ++z) if (!isdigit((unsigned char)rest[z])) { allDig = false; break; }
                if (allDig) dellen += std::atoi(rest.c_str() + up + 1);
            }
        }
        Variation& vref = getVariation(NIV, p, vn);
        int wustart = p - 200 > 1 ? p - 200 : 1;
        std::string wupseq = joinRef(ref, wustart, p - 1) + extra;
        int sanend = p + 200;
        std::string sanpseq = extra + joinRef(ref, p + dellen + (int)extra.size() - (int)extrains.size(), sanend);
        MismatchResult r3 = findMM3(vd, ref, p, sanpseq);
        MismatchResult r5 = findMM5(vd, ref, p + dellen + (int)extra.size() - (int)extrains.size() - 1, wupseq);

        std::vector<Mismatch> mmm = r3.mm; mmm.insert(mmm.end(), r5.mm.begin(), r5.mm.end());
        for (auto& mismatch : mmm) {
            std::string mm = mismatch.seq; int mp = mismatch.pos; int me = mismatch.end;
            if (mm.size() > 1) mm = std::string(1, mm[0]) + "&" + mm.substr(1);
            auto pit = NIV.find(mp); if (pit == NIV.end()) continue;
            auto vit = pit->second.find(mm); if (vit == pit->second.end()) continue;
            Variation& tv = vit->second;
            if (tv.varsCount == 0) continue;
            if (tv.meanQuality / tv.varsCount < cfg.goodq) continue;
            if (tv.meanPosition / tv.varsCount > (me == 3 ? r3.nm + 4 : r5.nm + 4)) continue;
            if (tv.varsCount >= dcnt + dellen || tv.varsCount / dcnt >= 8) continue;
            if (mp > p && me == 5) {
                double f = tv.meanPosition != 0 ? (mp - p) / (tv.meanPosition / (double)tv.varsCount) : 1;
                if (f > 1) f = 1;
                vd.refCoverage[p] += (int)(tv.varsCount * f);
                adjRefCnt(tv, ref.has(p) ? getVariationMaybe(NIV, p, ref.at(p)) : nullptr, dellen);
            }
            Variation* lref = (mp > p && me == 3 && ref.has(p)) ? getVariationMaybe(NIV, p, ref.at(p)) : nullptr;
            adjCnt(vref, tv, lref);
            pit->second.erase(mm);
            if (pit->second.empty()) NIV.erase(pit);
        }
        if (r3.misp != 0 && r3.mm.size() == 1) { auto pit = NIV.find(r3.misp);
            if (pit != NIV.end()) { auto vit = pit->second.find(r3.misnt);
                if (vit != pit->second.end() && vit->second.varsCount < dcnt) pit->second.erase(vit); } }
        if (r5.misp != 0 && r5.mm.size() == 1) { auto pit = NIV.find(r5.misp);
            if (pit != NIV.end()) { auto vit = pit->second.find(r5.misnt);
                if (vit != pit->second.end() && vit->second.varsCount < dcnt) pit->second.erase(vit); } }
        for (int sc5pp : r5.scp) {
            auto it = vd.softClips5End.find(sc5pp);
            if (it == vd.softClips5End.end() || it->second.used) continue;
            Sclip& tv = it->second;
            if (dcnt <= 2 && tv.varsCount / dcnt > 5) continue;
            std::string seq = findconseq(tv);
            if (!seq.empty() && ismatch(seq, wupseq, -1)) {
                if (sc5pp > p) vd.refCoverage[p] += tv.varsCount;
                adjCnt(vref, tv); tv.used = true;
            }
        }
        for (int sc3pp : r3.scp) {
            auto it = vd.softClips3End.find(sc3pp);
            if (it == vd.softClips3End.end() || it->second.used) continue;
            Sclip& tv = it->second;
            if (dcnt <= 2 && tv.varsCount / dcnt > 5) continue;
            std::string seq = findconseq(tv);
            if (!seq.empty() && ismatch(seq, substr(sanpseq, sc3pp - p), 1)) {
                if (sc3pp <= p) vd.refCoverage[p] += tv.varsCount;
                Variation* lref = (sc3pp <= p) ? nullptr : (ref.has(p) ? getVariationMaybe(NIV, p, ref.at(p)) : nullptr);
                adjCnt(vref, tv, lref); tv.used = true;
            }
        }
        // VariationRealigner.realigndel l.608-620 / vardict.pl l.4559: for a deletion with
        // microhomology where no read cleanly spans the gap as reference (noPassingReads), the
        // ambiguous reference reads at p (which only match the first bases of the repeat before
        // ending mid-gap) actually belong to the deletion, so fold h into vref.
        int pe = p + dellen + (int)extra.size() - (int)extrains.size();
        if (!bams.empty() && pe - p >= 5 && pe - p < maxReadLength - 10 && ref.has(p)) {
            Variation* h = getVariationMaybe(NIV, p, ref.at(p));
            if (h && h->varsCount != 0
                && noPassingReads(bams, region.chr, p, pe)
                && vref.varsCount > 2.0 * h->varsCount * (1 - (pe - p) / (double)maxReadLength)) {
                adjCnt(vref, *h, h);
            }
        }
    }
}

void realigndel(VariationData& vd, Reference& ref, const Config& cfg, const Region& region, int maxReadLength,
                const std::vector<BamReader*>& bams) {
    struct Item { int position; std::string desc; int count; };
    std::vector<Item> tmp;
    for (auto& [pos, m] : vd.positionToDeletionCount) for (auto& [d, c] : m) tmp.push_back({pos, d, c});
    std::sort(tmp.begin(), tmp.end(), [](const Item& a, const Item& b) {
        if (a.count != b.count) return a.count > b.count;
        if (a.position != b.position) return a.position < b.position;
        return a.desc > b.desc;
    });
    auto& NIV = vd.nonInsertionVariants;
    for (const auto& t : tmp) realignOneDel(vd, ref, cfg, region, maxReadLength, bams, t.position, t.desc, t.count);
    // Merge "-N&extra" into "-N".
    for (int i = (int)tmp.size() - 1; i > 0; --i) {
        int p = tmp[i].position; const std::string& vn = tmp[i].desc;
        auto pit = NIV.find(p); if (pit == NIV.end()) continue;
        auto vit = pit->second.find(vn); if (vit == pit->second.end()) continue;
        auto amp = vn.find('&');
        if (amp != std::string::npos && vn[0] == '-') {
            std::string tnkey = vn.substr(0, amp);
            auto tit = pit->second.find(tnkey);
            if (tit != pit->second.end() && vit->second.varsCount < tit->second.varsCount) {
                adjCnt(tit->second, vit->second); pit->second.erase(vit);
            }
        }
    }
}

// StructuralVariantsProcessor.adjSNV: merge a SHORT (<=5 bp) leftover soft-clip consensus into the
// adjacent SNV whose first base matches the clip (adjCnt sets extracnt + qstd, and coverage grows).
// This is the merge-back partner of CigarModifier's read-end soft-clipping; together they drive the
// ExtraAF / QStd columns.
void adjSNV(VariationData& vd, Reference& ref) {
    auto& NIV = vd.nonInsertionVariants;
    for (auto& [position, sclip] : vd.softClips5End) {
        if (sclip.used) continue;
        std::string seq = findconseq(sclip);
        if ((int)seq.size() > 5 || seq.empty()) continue;
        std::string bp = seq.substr(0, 1);
        int prev = position - 1;
        auto pit = NIV.find(prev);
        if (pit == NIV.end()) continue;
        auto vit = pit->second.find(bp);
        if (vit == pit->second.end()) continue;
        if (seq.size() > 1 && !(ref.has(position - 2) && ref.at(position - 2) == seq[1])) continue;
        adjCnt(vit->second, sclip);
        vd.refCoverage[prev] += sclip.varsCount;
    }
    for (auto& [position, sclip] : vd.softClips3End) {
        if (sclip.used) continue;
        std::string seq = findconseq(sclip);
        if ((int)seq.size() > 5 || seq.empty()) continue;
        std::string bp = seq.substr(0, 1);
        auto pit = NIV.find(position);
        if (pit == NIV.end()) continue;
        auto vit = pit->second.find(bp);
        if (vit == pit->second.end()) continue;
        if (seq.size() > 1 && !(ref.has(position + 1) && ref.at(position + 1) == seq[1])) continue;
        adjCnt(vit->second, sclip);
        vd.refCoverage[position] += sclip.varsCount;
    }
}

// ---- realignlgdel (large deletions from soft-clip breakpoints) -----------------------------------
// Faithful port of VariationRealigner.realignlgdel, including the bp==0 seed-based findMatch fallback
// and its partialPipeline coverage reload (via `reload`) when the breakpoint lands outside the region.
// The discordant-pair svcov/markSV bookkeeping is not ported (pairs=clusters=0; splits accumulates the
// clip count); this only affects the SV_info split/pair columns, not the emitted deletion's AF filter.

void realignlgdel(VariationData& vd, Reference& ref, const Config& cfg, const Region& region, int maxReadLength,
                  const SVReloadFn& reload, const std::vector<BamReader*>& bams) {
    auto& NIV = vd.nonInsertionVariants;
    const int EXT = Config::EXTENSION;

    auto collectSorted = [&](std::map<int, Sclip>& clips) {
        std::vector<std::pair<int, Sclip*>> v;
        for (auto& [p, sc] : clips)
            if (p >= region.start - EXT && p <= region.end + EXT) v.push_back({p, &sc});
        std::sort(v.begin(), v.end(), [](auto& a, auto& b) {
            if (a.second->varsCount != b.second->varsCount) return a.second->varsCount > b.second->varsCount;
            return a.first < b.first;
        });
        return v;
    };

    // 5' soft-clipped reads
    for (auto& [p, sc5vp] : collectSorted(vd.softClips5End)) {
        Sclip& sc5v = *sc5vp;
        int cnt = sc5v.varsCount;
        if (cnt < cfg.minReads) break;
        if (sc5v.used) continue;
        std::string seq = findconseq(sc5v);
        if (seq.empty() || (int)seq.size() < 7) continue;
        int bp = findbp(seq, p - 5, ref, -1, vd.chrLen, cfg.indelsize);
        std::string EXTRA;
        if (bp == 0) {   // findbp failed: try seed-based findMatch (SV-cluster / partialPipeline gated)
            if (islowcomplexseq(seq)) continue;
            Match match = findMatch(seq, ref, p, -1, Reference::SEED_1, 1);
            bp = match.bp; EXTRA = match.extra;
            if (!(bp != 0 && p - bp > 15 && p - bp < Config::SVMAXLEN)) continue;
            bp++;
            if (cnt <= cfg.minReads) continue; // svcov==0 path (no discordant-pair support)
            // SV split-read marker (VariationRealigner ~1042-1046). markSV (discordant-pair
            // clusters) is not ported -> pairs=clusters=0; splits accumulates the clip count.
            { SVInfo& sv = vd.svInfoAt[bp]; sv.type = "DEL"; sv.splits += cnt; }
            // partialPipeline reload (VariationRealigner 1048-1061): if the breakpoint lands before
            // the region, re-read coverage at [bp-maxReadLength, min(bp+maxReadLength, region.start-1)]
            // so refCoverage[bp] reflects the true depth and a low-VAF deletion is AF-filtered.
            if (bp < region.start) {
                int tts = bp - maxReadLength;
                int tte = (bp + maxReadLength >= region.start) ? region.start - 1 : bp + maxReadLength;
                reload(tts + 200, tte - 200);
            }
        }
        int dellen = p - bp;
        std::string extra;
        std::string gt = "-" + std::to_string(dellen);
        if (EXTRA.empty()) {
            int en = 0;
            while (en < (int)seq.size() && !(ref.has(bp - en - 1) && seq[en] == ref.at(bp - en - 1))) {
                extra += seq[en]; en++;
            }
            if (!extra.empty()) { extra = reverseStr(extra); gt = "-" + std::to_string(dellen) + "&" + extra; bp -= (int)extra.size(); }
        } else {
            dellen -= (int)EXTRA.size();
            gt = dellen == 0 ? "-" + std::to_string((int)EXTRA.size()) + "^" + EXTRA
                             : "-" + std::to_string(dellen) + "&" + EXTRA;
        }
        // 3' breakpoint (sc3p)
        int n = 0;
        if (extra.empty() && EXTRA.empty()) {
            while (ref.has(bp + n) && ref.has(bp + dellen + n) && ref.at(bp + n) == ref.at(bp + dellen + n)) n++;
        }
        int sc3p = bp + n;
        std::string str; int mcnt = 0;
        while (mcnt <= 3 && ref.has(bp + n) && ref.has(bp + dellen + n) && ref.at(bp + n) != ref.at(bp + dellen + n)) {
            str += ref.at(bp + dellen + n); n++; mcnt++;
        }
        if (str.size() == 1) {
            int nm = 0;
            while (ref.has(bp + n) && ref.has(bp + dellen + n) && ref.at(bp + n) == ref.at(bp + dellen + n)) { n++; if (n != 0) nm++; }
            if (nm >= 3 && !vd.softClips3End.count(sc3p)) sc3p = bp + n;
        }
        Variation& tv = getVariation(NIV, bp, gt);
        tv.qstd = true; tv.pstd = true;
        adjCnt(tv, sc5v); sc5v.used = (bp != 0);
        if (!vd.refCoverage.count(bp) && vd.refCoverage.count(p)) vd.refCoverage[bp] = vd.refCoverage[p];
        if (dellen < cfg.indelsize) for (int tp = bp; tp < bp + dellen; ++tp) vd.refCoverage[tp] += sc5v.varsCount;
        auto s3it = vd.softClips3End.find(sc3p);
        if (s3it != vd.softClips3End.end() && !s3it->second.used) {
            Sclip& sclip = s3it->second;
            if (sc3p > bp) adjCnt(tv, sclip, ref.has(bp) ? getVariationMaybe(NIV, bp, ref.at(bp)) : nullptr);
            else adjCnt(tv, sclip);
            if (sc3p == bp && dellen < cfg.indelsize)
                for (int tp = bp; tp < bp + dellen; ++tp) vd.refCoverage[tp] += sclip.varsCount;
            for (int ip = bp + 1; ip < sc3p; ++ip) {
                if (!ref.has(dellen + ip)) continue;
                std::string rk(1, ref.at(dellen + ip));
                auto pit = NIV.find(ip);
                if (pit == NIV.end()) continue;
                auto vit = pit->second.find(rk);
                if (vit == pit->second.end()) continue;
                rmCnt(vit->second, sclip);
                if (vit->second.varsCount == 0) pit->second.erase(vit);
                if (pit->second.empty()) NIV.erase(pit);
            }
            sclip.used = (bp != 0);
        }
        // Java realignlgdel l.1173-1179: re-run realigndel on this single large deletion (mismatch-SNV
        // attribution + noPassingReads fold-in) and bump the SV anchor's split count by reads it drew in.
        {
            int origCnt = tv.varsCount;
            realignOneDel(vd, ref, cfg, region, maxReadLength, bams, bp, gt, origCnt);
            auto pit2 = NIV.find(bp);
            int newCnt = (pit2 != NIV.end() && pit2->second.count(gt)) ? pit2->second[gt].varsCount : origCnt;
            auto svit = vd.svInfoAt.find(bp);
            if (svit != vd.svInfoAt.end()) svit->second.splits += newCnt - origCnt;
        }
    }

    // 3' soft-clipped reads
    for (auto& [p, sc3vp] : collectSorted(vd.softClips3End)) {
        Sclip& sc3v = *sc3vp;
        int cnt = sc3v.varsCount;
        if (cnt < cfg.minReads) break;
        if (sc3v.used) continue;
        std::string seq = findconseq(sc3v);
        if (seq.empty() || (int)seq.size() < 7) continue;
        int bp = findbp(seq, p + 5, ref, 1, vd.chrLen, cfg.indelsize);
        std::string EXTRA;
        if (bp == 0) {
            if (islowcomplexseq(seq)) continue;
            Match match = findMatch(seq, ref, p, 1, Reference::SEED_1, 1);
            bp = match.bp; EXTRA = match.extra;
            if (!(bp != 0 && bp - p > 15)) continue;
            if (cnt <= cfg.minReads) continue;
            // SV split-read marker (VariationRealigner ~1254-1258). markSV not ported ->
            // pairs=clusters=0; marker is keyed at p (the 3' clip position == variant bp).
            { SVInfo& sv = vd.svInfoAt[p]; sv.type = "DEL"; sv.splits += cnt; }
            // partialPipeline reload (VariationRealigner 1260-1273): breakpoint past the region end ->
            // re-read coverage at [max(bp-maxReadLength, region.end+1), bp+maxReadLength].
            if (bp > region.end) {
                int tts = (bp - maxReadLength <= region.end) ? region.end + 1 : bp - maxReadLength;
                int tte = bp + maxReadLength;
                reload(tts + 200, tte - 200);
            }
        }
        int dellen = bp - p;
        std::string extra;
        if (!EXTRA.empty()) dellen -= (int)EXTRA.size();
        else { int en = 0; while (en < (int)seq.size() && !(ref.has(bp + en) && seq[en] == ref.at(bp + en))) { extra += seq[en]; en++; } }
        std::string gt = "-" + std::to_string(dellen);
        bp = p; // set to 5'
        if (!extra.empty()) gt = "-" + std::to_string(dellen) + "&" + extra;
        else if (!EXTRA.empty()) gt = "-" + std::to_string(dellen) + "&" + EXTRA;
        else {
            while (ref.has(bp - 1) && ref.has(bp + dellen - 1) && ref.at(bp - 1) == ref.at(bp + dellen - 1)) bp--;
            // The 5' walk moved the breakpoint off p; carry the SV split marker with it
            // (VariationRealigner ~1303-1310).
            if (bp != p) {
                auto svit = vd.svInfoAt.find(p);
                if (svit != vd.svInfoAt.end()) { vd.svInfoAt[bp] = svit->second; vd.svInfoAt.erase(svit); }
            }
        }
        Variation& tv = getVariation(NIV, bp, gt);
        tv.qstd = true; tv.pstd = true;
        if (dellen < cfg.indelsize)
            for (int tp = bp; tp < bp + dellen + (int)extra.size(); ++tp) vd.refCoverage[tp] += sc3v.varsCount;
        if (!vd.refCoverage.count(bp)) vd.refCoverage[bp] = vd.refCoverage.count(p - 1) ? vd.refCoverage[p - 1] : sc3v.varsCount;
        sc3v.meanPosition += (double)dellen * sc3v.varsCount;
        adjCnt(tv, sc3v); sc3v.used = true;
        // Java realignlgdel l.1341-1348: re-run realigndel on this single large deletion.
        {
            int origCnt = tv.varsCount;
            realignOneDel(vd, ref, cfg, region, maxReadLength, bams, bp, gt, origCnt);
            auto pit2 = NIV.find(bp);
            int newCnt = (pit2 != NIV.end() && pit2->second.count(gt)) ? pit2->second[gt].varsCount : origCnt;
            auto svit = vd.svInfoAt.find(bp);
            if (svit != vd.svInfoAt.end()) svit->second.splits += newCnt - origCnt;
        }
    }
}

// ---- realignlgins30 (large insertions from paired 5'/3' soft-clips) ------------------------------

struct Match35 { int b5; int b3; int score; };
static Match35 find35match(const std::string& seq5, const std::string& seq3) {
    const int longMismatch = 2;
    int maxLen = 0, b3 = 0, b5 = 0;
    const int s5 = (int)seq5.size(), s3 = (int)seq3.size();
    for (int i = 0; i < s5 - 8; ++i) {
        for (int j = 1; j < s3 - 8; ++j) {
            int nm = 0, n = 0;
            while (n + j <= s3 && i + n <= s5) {
                // substr(seq3,-j-n,1) is the char at s3-(j+n) (always in range here); substr(seq5,i+n,1)
                // is seq5[i+n], or "" when i+n==s5 -> represented by a sentinel that always mismatches.
                char c3 = seq3[s3 - (j + n)];
                char c5 = (i + n < s5) ? seq5[i + n] : (char)0xFF;
                if (c3 != c5) nm++;
                if (nm > longMismatch) break;
                n++;
            }
            if (n - nm > maxLen && n - nm > 8 && nm / (double)n < 0.1 &&
                (n + j >= (int)seq3.size() || i + n >= (int)seq5.size())) {
                return { i, j, n - nm };
            }
        }
    }
    return { b5, b3, maxLen };
}

void realignlgins30(VariationData& vd, Reference& ref, const Config& cfg, const Region& region, int maxReadLength) {
    auto& NIV = vd.nonInsertionVariants;
    const int EXT = Config::EXTENSION;
    auto collect = [&](std::map<int, Sclip>& clips) {
        std::vector<std::pair<int, Sclip*>> v;
        for (auto& [p, sc] : clips)
            if (p >= region.start - EXT && p <= region.end + EXT) v.push_back({p, &sc});
        std::sort(v.begin(), v.end(), [](auto& a, auto& b) {
            if (a.second->varsCount != b.second->varsCount) return a.second->varsCount > b.second->varsCount;
            return a.first < b.first;
        });
        return v;
    };
    auto tmp5 = collect(vd.softClips5End);
    auto tmp3 = collect(vd.softClips3End);
    for (auto& [p5, sc5vp] : tmp5) {
        Sclip& sc5v = *sc5vp;
        int cnt5 = sc5v.varsCount;
        if (sc5v.used) continue;
        for (auto& [p3, sc3vp] : tmp3) {
            Sclip& sc3v = *sc3vp;
            int cnt3 = sc3v.varsCount;
            if (sc5v.used) break;
            if (sc3v.used) continue;
            if (p5 - p3 > maxReadLength * 2.5) continue;
            if (p3 - p5 > maxReadLength - 10) continue;
            std::string seq5 = findconseq(sc5v), seq3 = findconseq(sc3v);
            if ((int)seq5.size() <= 10 || (int)seq3.size() <= 10) continue;
            if (!(cnt5 / (double)cnt3 >= 0.08 && cnt5 / (double)cnt3 <= 12)) continue;
            Match35 m = find35match(seq5, seq3);
            int bp5 = m.b5, bp3 = m.b3, score = m.score;
            if (score == 0) continue;
            int smscore = score / 2;
            std::string ins = (bp3 + smscore > 1) ? substr(seq3, 0, -(bp3 + smscore) + 1) : seq3;
            if (bp5 + smscore > 0) ins += reverseStr(substr(seq5, 0, bp5 + smscore));
            if (islowcomplexseq(ins)) continue;
            int bi = 0; bool isIns = false, isDel = false; std::string key;
            if (p5 > p3) {
                if ((int)seq3.size() > (int)ins.size() &&
                    !ismatch(substr(seq3, (int)ins.size()), joinRef(ref, p5, p5 + (int)seq3.size() - (int)ins.size() + 2), 1)) continue;
                if ((int)seq5.size() > (int)ins.size() &&
                    !ismatch(substr(seq5, (int)ins.size()), joinRef(ref, p3 - ((int)seq5.size() - (int)ins.size()) - 2, p3 - 1), -1)) continue;
                std::string tmp = joinRef(ref, p3, p5 - 1);
                if ((int)tmp.size() > (int)ins.size()) { key = std::to_string(p3 - p5) + "^" + ins; bi = p3; isDel = true; }
                else if ((int)tmp.size() < (int)ins.size()) {
                    ins = substr(ins, 0, (int)ins.size() - (int)tmp.size()) + "&" + substr(ins, p3 - p5);
                    key = "+" + ins; bi = p3 - 1; isIns = true;
                } else { key = "-" + std::to_string((int)ins.size()) + "^" + ins; bi = p3; isDel = true; }
            } else {
                if ((int)seq3.size() > (int)ins.size() &&
                    !ismatch(substr(seq3, (int)ins.size()), joinRef(ref, p5, p5 + (int)seq3.size() - (int)ins.size() + 2), 1)) continue;
                if ((int)seq5.size() > (int)ins.size() &&
                    !ismatch(substr(seq5, (int)ins.size()), joinRef(ref, p3 - ((int)seq5.size() - (int)ins.size()) - 2, p3 - 1), -1)) continue;
                std::string tmp;
                if ((int)ins.size() <= p3 - p5) {
                    int rpt = 2, tnr = 3;
                    while (((p3 - p5 + (int)ins.size()) / (double)tnr) / (double)ins.size() > 1) {
                        if ((p3 - p5 + (int)ins.size()) % tnr == 0) rpt++;
                        tnr++;
                    }
                    tmp += joinRef(ref, p5, p5 + (p3 - p5 + (int)ins.size()) / (double)rpt - (int)ins.size());
                    ins = "+" + tmp + ins;
                } else {
                    tmp += joinRef(ref, p5, p3 - 1);
                    if (((int)ins.size() - (int)tmp.size()) % 2 == 0) {
                        int tex = ((int)ins.size() - (int)tmp.size()) / 2;
                        ins = (tmp + substr(ins, 0, tex)) == substr(ins, tex) ? ("+" + substr(ins, tex)) : "+" + tmp + ins;
                    } else ins = "+" + tmp + ins;
                }
                key = ins; bi = p5 - 1; isIns = true;
            }
            Variation& vref = isIns ? vd.insertionVariants[bi][key] : NIV[bi][key];
            sc3v.used = true; sc5v.used = true;
            vref.pstd = true; vref.qstd = true;
            vd.refCoverage[bi] += sc5v.varsCount;
            if (isIns) {
                Variation* mvref = ref.has(bi) ? getVariationMaybe(NIV, bi, ref.at(bi)) : nullptr;
                adjCnt(vref, sc3v, mvref); adjCnt(vref, sc5v);
                // VariationRealigner.realignlgins30 re-runs realignins on the allele it just created
                // (Java: realignins({bi:{ins:vref.varsCount}})); the whole-map realignins pass ran
                // earlier in the pipeline and never sees this insertion. This attracts the neighbouring
                // mismatched reads that make up the insertion's per-variant mean-mismatch (NM) count.
                realignOneIns(vd, ref, cfg, maxReadLength, bi, key, vref.varsCount);
            } else if (isDel) {
                adjCnt(vref, sc3v, ref.has(bi) ? getVariationMaybe(NIV, bi, ref.at(bi)) : nullptr);
                adjCnt(vref, sc5v);
            } else { adjCnt(vref, sc3v); adjCnt(vref, sc5v); }
            break;
        }
    }
}

// ---- realignlgins (tandem insertions / duplications from soft-clips) -----------------------------

struct BaseInsertion { int bi; std::string ins; int bi2; };
static BaseInsertion adjInsPosR(int bi, std::string ins, Reference& ref) {
    int n = 1, len = (int)ins.size();
    while (ref.has(bi) && ref.at(bi) == ins[ins.size() - n]) { n++; if (n > len) n = 1; bi--; }
    if (n > 1) ins = substr(ins, 1 - n) + substr(ins, 0, 1 - n);
    return { bi, ins, bi };
}

// VariationRealigner.findbi: slide the soft-clip consensus against the reference to detect a tandem
// insertion (a shift where the overhang re-matches). Returns {breakpoint, inserted seq, bi2}.
static BaseInsertion findbi(const std::string& seq, int position, Reference& ref, int dir, int chrLen) {
    const int maxmm = 3;
    int dirExt = dir == -1 ? 1 : 0;
    int score = 0, bi = 0, bi2 = 0;
    std::string ins;
    for (int n = 6; n < (int)seq.size(); ++n) {
        if (position + 6 >= chrLen) break;
        int mm = 0, i = 0;
        std::set<char> m;
        for (i = 0; i + n < (int)seq.size(); ++i) {
            int rp = position + dir * i - dirExt;
            if (rp < 1 || rp > chrLen) break;
            if (!(ref.has(rp) && seq[i + n] == ref.at(rp))) mm++;
            else m.insert(seq[i + n]);
            if (mm > maxmm) break;
        }
        int mnt = (int)m.size();
        if (mnt < 2) continue;
        if ((mnt >= 3 && i + n >= (int)seq.size() - 1 && i >= 8 && mm / (double)i < 0.15) ||
            (mnt >= 2 && mm == 0 && i + n == (int)seq.size() && n >= 20 && i >= 8)) {
            std::string insert = substr(seq, 0, n), extra;
            int ept = 0;
            while (n + ept + 1 < (int)seq.size() &&
                   (!(ref.has(position + ept * dir - dirExt) && seq[n + ept] == ref.at(position + ept * dir - dirExt)) ||
                    !(ref.has(position + (ept + 1) * dir - dirExt) && seq[n + ept + 1] == ref.at(position + (ept + 1) * dir - dirExt)))) {
                extra += seq[n + ept]; ept++;
            }
            if (dir == -1) {
                insert += extra;
                std::reverse(insert.begin(), insert.end());
                if (!extra.empty()) insert.insert(insert.size() - extra.size(), "&");
                if (mm == 0 && i + n == (int)seq.size()) {
                    bi = position - 1 - (int)extra.size(); ins = insert; bi2 = position - 1;
                    if (extra.empty()) return adjInsPosR(bi, ins, ref);
                    return { bi, ins, bi2 };
                } else if (i - mm > score) { bi = position - 1 - (int)extra.size(); ins = insert; bi2 = position - 1; score = i - mm; }
            } else {
                int s = -1;
                if (!extra.empty()) insert += "&" + extra;
                else {
                    while (s >= -n && charAt(insert, s) != (char)-1 && ref.has(position + s) && charAt(insert, s) == ref.at(position + s)) s--;
                    if (s < -1) {
                        std::string tins = substr(insert, s + 1, 1 - s);
                        insert.erase(insert.size() + s + 1);
                        insert = tins + insert;
                    }
                }
                if (mm == 0 && i + n == (int)seq.size()) {
                    bi = position + s; ins = insert; bi2 = position + s + (int)extra.size();
                    if (extra.empty()) return adjInsPosR(bi, ins, ref);
                    return { bi, ins, bi2 };
                } else if (i - mm > score) { bi = position + s; ins = insert; bi2 = position + s + (int)extra.size(); score = i - mm; }
            }
        }
    }
    if (bi2 == bi && !ins.empty() && bi != 0) return adjInsPosR(bi, ins, ref);
    return { bi, ins, bi2 };
}

// StructuralVariantsProcessor.isOverlap: shared by markSVDel/markDUPSV. (Definition relocated here so
// the split-read DUP path in realignlgins can reuse it; markSVDel below still calls it.)
static bool isOverlapSV(int s1, int e1, int s2, int e2, int rlen) {
    if (s1 >= e2 || s2 >= e1) return false;
    int p[4] = {s1, e1, s2, e2}; std::sort(p, p + 4);
    int ins = p[2] - p[1];
    if (e1 != s1 && e2 != s2 && ins / (double)(e1 - s1) > 0.75 && ins / (double)(e2 - s2) > 0.75) return true;
    if ((p[1] - p[0]) + (p[3] - p[2]) < 3 * rlen) return true;
    return false;
}

// StructuralVariantsProcessor.markDUPSV: mark same-orientation dup clusters overlapping [start,end]
// used; returns {count of clusters marked (clusters), sum of their varsCount (pairs)}. Note the
// start2/end2 convention differs from markSV (start/mend vs mstart/end swapped).
static std::pair<int,int> markDUPSV(int start, int end,
                                    std::initializer_list<std::vector<Sclip>*> lists, int rlen) {
    int pairs = 0, cnt = 0;
    for (auto* lst : lists) {
        for (auto& sv_r : *lst) {
            int s2, e2;
            if (sv_r.start < sv_r.mstart) { s2 = sv_r.start; e2 = sv_r.mend; }
            else                          { s2 = sv_r.mstart; e2 = sv_r.end; }
            if (isOverlapSV(start, end, s2, e2, rlen)) {
                sv_r.used = true;
                cnt++;
                pairs += sv_r.varsCount;
            }
        }
    }
    return {cnt, pairs};
}

void realignlgins(VariationData& vd, Reference& ref, const Config& cfg, const Region& region,
                  int maxReadLength, const SVReloadFn& reload) {
    auto& NIV = vd.nonInsertionVariants;
    const int EXT = Config::EXTENSION;
    auto collect = [&](std::map<int, Sclip>& clips) {
        std::vector<std::pair<int, Sclip*>> v;
        for (auto& [p, sc] : clips) if (p >= region.start - EXT && p <= region.end + EXT) v.push_back({p, &sc});
        std::sort(v.begin(), v.end(), [](auto& a, auto& b) {
            if (a.second->varsCount != b.second->varsCount) return a.second->varsCount > b.second->varsCount;
            return a.first < b.first;
        });
        return v;
    };
    // 5' soft-clips
    for (auto& [p, sc5vp] : collect(vd.softClips5End)) {
        Sclip& sc5v = *sc5vp;
        if (sc5v.varsCount < cfg.minReads) break;
        if (sc5v.used) continue;
        std::string seq = findconseq(sc5v);
        if (seq.empty() || (int)seq.size() < 12) continue;
        BaseInsertion tpl = findbi(seq, p, ref, -1, vd.chrLen);
        int bi = tpl.bi; std::string ins = tpl.ins;
        bool madeSV = false;
        if (bi == 0) {   // findbi failed: seed findMatch DUP path (markDUPSV / partialPipeline gated)
            if (islowcomplexseq(seq)) continue;
            Match match = findMatch(seq, ref, p, -1, Reference::SEED_1, 1);
            bi = match.bp; std::string EXTRA = match.extra;
            if (!(bi != 0 && bi - p > 15 && bi - p < Config::SVMAXLEN)) continue;
            // Large insertion whose partner lands past the region end: re-read coverage over
            // [tts, tte] so refCoverage reflects the true depth (VariationRealigner 1648-1662).
            if (bi > region.end) {
                int tts = bi - maxReadLength;
                int tte = bi + maxReadLength;
                if (bi - maxReadLength <= region.end) tts = region.end + 1;
                reload(tts + 200, tte - 200);
            }
            if (bi - p > cfg.SVMINLEN + 2 * Config::SVFLANK) {
                ins = joinRef(ref, p, p + Config::SVFLANK - 1);
                ins += "<dup" + std::to_string(bi - p - 2 * Config::SVFLANK + 1) + ">";
                ins += joinRefFor5Lgins(ref, bi - Config::SVFLANK + 1, bi, seq, EXTRA);
            } else {
                ins = joinRefFor5Lgins(ref, p, bi, seq, EXTRA);
            }
            ins += EXTRA;
            // markDUPSV(p, bi): mark same-orientation dup clusters overlapping [p, bi] and collect
            // their pair support (VariationRealigner 1673). bi here is still the match position.
            auto dupmark = markDUPSV(p, bi, {&vd.svfdup, &vd.svrdup}, maxReadLength);
            if (!vd.refCoverage.count(p - 1) || (vd.refCoverage.count(bi) && vd.refCoverage[p - 1] < vd.refCoverage[bi])) {
                vd.refCoverage[p - 1] = vd.refCoverage.count(bi) ? vd.refCoverage[bi] : sc5v.varsCount;
            } else if (sc5v.varsCount > vd.refCoverage[p - 1]) vd.refCoverage[p - 1] += sc5v.varsCount;
            bi = p - 1;
            // VariationRealigner.getSV(nonInsertionVariants, bi): anchor a DUP SV marker at bi so
            // ToVarsBuilder (which iterates nonInsertionVariants) visits this position and emits the
            // insertion. splits accumulates the split count; pairs/clusters come from markDUPSV.
            NIV[bi];  // create the non-insertion anchor (empty allele map)
            SVInfo& sv = vd.svInfoAt[bi];
            sv.type = "DUP"; sv.splits += sc5v.varsCount;
            sv.pairs += dupmark.second; sv.clusters += dupmark.first;
            madeSV = true;
        }
        Variation& iref = vd.insertionVariants[bi]["+" + ins];
        iref.pstd = true; iref.qstd = true;
        adjCnt(iref, sc5v);
        // Java increments refCoverage only when the anchor carries no SV marker (sv == null).
        if (NIV.count(bi) && !madeSV && !vd.svInfoAt.count(bi)) vd.refCoverage[bi] += sc5v.varsCount;
        int len = (int)ins.size(); if (ins.find('&') != std::string::npos) len--;
        if (!sc5v.seq.empty()) {
            int seqLen = sc5v.seq.rbegin()->first + 1;
            for (int ii = len + 1; ii < seqLen; ++ii) {
                int pii = bi - ii + len;
                auto sit = sc5v.seq.find(ii); if (sit == sc5v.seq.end()) continue;
                for (auto& [tnt, tv] : sit->second) {
                    Variation& tvr = getVariation(NIV, pii, std::string(1, tnt));
                    adjCnt(tvr, *tv); tvr.pstd = true; tvr.qstd = true;
                    vd.refCoverage[pii] += tv->varsCount;
                }
            }
        }
        sc5v.used = (bi + len != 0);
        // VariationRealigner.realignlgins (1738-1748): re-run realignins on the freshly created
        // insertion (the whole-map pass ran before this insertion existed) so it attracts the reads
        // that span it, then bump the SV anchor's split count by the newly attracted reads.
        {
            int origCount = iref.varsCount;
            realignOneIns(vd, ref, cfg, maxReadLength, bi, "+" + ins, origCount);
            Variation& kref = vd.insertionVariants[bi]["+" + ins];
            auto svit = vd.svInfoAt.find(bi);
            if (svit != vd.svInfoAt.end()) svit->second.splits += kref.varsCount - origCount;
        }
    }
    // 3' soft-clips
    for (auto& [p, sc3vp] : collect(vd.softClips3End)) {
        Sclip& sc3v = *sc3vp;
        if (sc3v.varsCount < cfg.minReads) break;
        if (sc3v.used) continue;
        std::string seq = findconseq(sc3v);
        if (seq.empty() || (int)seq.size() < 12) continue;
        BaseInsertion tpl = findbi(seq, p, ref, 1, vd.chrLen);
        int bi = tpl.bi; std::string ins = tpl.ins;
        if (bi == 0) {
            if (islowcomplexseq(seq)) continue;
            Match match = findMatch(seq, ref, p, 1, Reference::SEED_1, 1);
            bi = match.bp; std::string EXTRA = match.extra;
            if (!(bi != 0 && p - bi > 15 && p - bi < Config::SVMAXLEN)) continue;
            // Large insertion whose partner lands before the region start: re-read coverage over
            // [tts, tte] so refCoverage reflects the true depth (VariationRealigner 1812-1828).
            if (bi < region.start) {
                int tts = bi - maxReadLength;
                int tte = bi + maxReadLength;
                if (bi + maxReadLength >= region.start) tte = region.start - 1;
                reload(tts + 200, tte - 200);
            }
            int shift5 = 0;
            while (ref.has(p - 1) && ref.has(bi - 1) && ref.at(p - 1) == ref.at(bi - 1)) { p--; bi--; shift5++; }
            if (p - bi > cfg.SVMINLEN + 2 * Config::SVFLANK) {
                ins = joinRefFor3Lgins(ref, bi, bi + Config::SVFLANK - 1, shift5, seq, EXTRA);
                ins += "<dup" + std::to_string(p - bi - 2 * Config::SVFLANK) + ">";
                ins += joinRef(ref, p - Config::SVFLANK, p - 1);
            } else {
                ins = joinRefFor3Lgins(ref, bi, p - 1, shift5, seq, EXTRA);
            }
            ins += EXTRA;
            // markDUPSV(bi, p-1): same-orientation dup cluster pair support (VariationRealigner 1847).
            auto dupmark = markDUPSV(bi, p - 1, {&vd.svfdup, &vd.svrdup}, maxReadLength);
            bi = bi - 1;
            // getSV anchor (see 5' branch): make ToVarsBuilder visit bi and emit the insertion.
            NIV[bi];
            SVInfo& sv = vd.svInfoAt[bi];
            sv.type = "DUP"; sv.splits += sc3v.varsCount;
            sv.pairs += dupmark.second; sv.clusters += dupmark.first;
            if (!vd.refCoverage.count(bi) || (vd.refCoverage.count(p) && vd.refCoverage[bi] < vd.refCoverage[p])) {
                vd.refCoverage[bi] = vd.refCoverage.count(p) ? vd.refCoverage[p] : sc3v.varsCount;
            } else if (sc3v.varsCount > vd.refCoverage[bi]) vd.refCoverage[bi] += sc3v.varsCount;
        }
        Variation& iref = vd.insertionVariants[bi]["+" + ins];
        iref.pstd = true; iref.qstd = true;
        Variation* lref = ref.has(bi) ? getVariationMaybe(NIV, bi, ref.at(bi)) : nullptr;
        double m3 = sc3v.varsCount ? sc3v.meanPosition / sc3v.varsCount : 0;
        if (p - bi > m3) lref = nullptr;
        adjCnt(iref, sc3v, lref);
        int len = (int)ins.size(); if (ins.find('&') != std::string::npos) len--;
        if (!sc3v.seq.empty()) {
            int lenSeq = sc3v.seq.rbegin()->first + 1;
            for (int ii = len; ii < lenSeq; ++ii) {
                int pii = p + ii - len;
                auto sit = sc3v.seq.find(ii); if (sit == sc3v.seq.end()) continue;
                for (auto& [tnt, tv] : sit->second) {
                    Variation& vref = getVariation(NIV, pii, std::string(1, tnt));
                    adjCnt(vref, *tv); vref.pstd = true; vref.qstd = true;
                    vd.refCoverage[pii] += tv->varsCount;
                }
            }
        }
        sc3v.used = true;
        // VariationRealigner.realignlgins (1917-1920): re-run realignins on the created insertion and
        // bump the SV anchor's split count by the newly attracted reads (see the 5' branch).
        {
            int origCount = iref.varsCount;
            realignOneIns(vd, ref, cfg, maxReadLength, bi, "+" + ins, origCount);
            Variation& kref = vd.insertionVariants[bi]["+" + ins];
            auto svit = vd.svInfoAt.find(bi);
            if (svit != vd.svInfoAt.end()) svit->second.splits += kref.varsCount - origCount;
        }
    }
}

// ---- Structural variants: discordant-pair deletions ---------------------------------------------
// Ports VariationRealigner.filterSV/checkCluster + StructuralVariantsProcessor.findDELdisc/markSV/
// isOverlap for the DEL discordant-pair path. Turns svfdel/svrdel clusters into a <DEL> Variation at
// the breakpoint plus an SVInfo (splits/pairs/clusters) marker consumed by callVariants.

namespace {
struct Cluster {
    int mateStart_ms=0, mateEnd_me=0, cnt=0, mateLength_mlen=0, start_s=0, end_e=0;
    double pmean_rp=0, qmean_q=0, Qmean_Q=0, nm=0;
};

// VariationRealigner.checkCluster: collapse a cluster's mates into the dominant sub-cluster (mates
// grouped by mate start within MINSVCDIST*rlen). Returns an all-zero Cluster (mateStart_ms==0) when
// the top sub-cluster holds < 60% of the mates.
Cluster checkCluster(std::vector<Mate>& mates, int rlen) {
    std::stable_sort(mates.begin(), mates.end(),
                     [](const Mate& a, const Mate& b) { return a.mateStart_ms < b.mateStart_ms; });
    std::vector<Cluster> clusters;
    { const Mate& f = mates[0]; Cluster c; c.mateStart_ms=f.mateStart_ms; c.mateEnd_me=f.mateEnd_me;
      c.start_s=f.start_s; c.end_e=f.end_e; clusters.push_back(c); }
    int cur = 0;
    for (const Mate& m : mates) {
        if (m.mateStart_ms - clusters[cur].mateEnd_me > Config::MINSVCDIST * rlen) {
            Cluster c; c.mateStart_ms=m.mateStart_ms; c.mateEnd_me=m.mateEnd_me;
            c.start_s=m.start_s; c.end_e=m.end_e; clusters.push_back(c); cur++;
        }
        Cluster& g = clusters[cur];
        g.cnt++;
        g.mateLength_mlen += m.mateLength_mlen;
        if (m.mateEnd_me > g.mateEnd_me) g.mateEnd_me = m.mateEnd_me;
        if (m.start_s < g.start_s) g.start_s = m.start_s;
        if (m.end_e > g.end_e) g.end_e = m.end_e;
        g.pmean_rp += m.pmean_rp; g.qmean_q += m.qmean_q; g.Qmean_Q += m.Qmean_Q; g.nm += m.nm;
    }
    std::stable_sort(clusters.begin(), clusters.end(),
                     [](const Cluster& a, const Cluster& b) { return a.cnt > b.cnt; });
    const Cluster& fc = clusters[0];
    Cluster res;
    if (fc.cnt / (double)mates.size() >= 0.60) {
        res = fc;
        res.mateLength_mlen = fc.mateLength_mlen / fc.cnt;
    }
    return res;
}

void filterSVList(std::vector<Sclip>& list, int maxReadLength) {
    for (auto& sv : list) {
        if (sv.mates.empty()) { sv.used = true; continue; }
        Cluster cl = checkCluster(sv.mates, maxReadLength);
        if (cl.mateStart_ms != 0) {
            sv.mstart = cl.mateStart_ms; sv.mend = cl.mateEnd_me; sv.varsCount = cl.cnt;
            sv.mlen = cl.mateLength_mlen; sv.start = cl.start_s; sv.end = cl.end_e;
            sv.meanPosition = cl.pmean_rp; sv.meanQuality = cl.qmean_q;
            sv.meanMappingQuality = cl.Qmean_Q; sv.numberOfMismatches = cl.nm;
        } else {
            sv.used = true;
        }
        if (sv.disc != 0 && sv.varsCount / (double)sv.disc < 0.5) {
            if (!(sv.varsCount / (double)sv.disc >= 0.35 && sv.varsCount >= 5)) sv.used = true;
        }
        // dominant soft-clip position (max count; lowest position on tie)
        int bestp = 0, bestc = -1;
        for (auto& [p, cnt] : sv.soft) if (cnt > bestc) { bestc = cnt; bestp = p; }
        sv.softp = sv.soft.empty() ? 0 : bestp;
    }
}

// StructuralVariantsProcessor.markSV: mark reciprocal-orientation clusters overlapping [start,end] used.
void markSVDel(int start, int end, std::vector<Sclip>& list, int rlen) {
    for (auto& sv_r : list) {
        int s2, e2;
        if (sv_r.start < sv_r.mstart) { s2 = sv_r.end; e2 = sv_r.mstart; }
        else                          { s2 = sv_r.mend; e2 = sv_r.start; }
        if (isOverlapSV(start, end, s2, e2, rlen)) sv_r.used = true;
    }
}
} // anonymous namespace

void filterSVStructures(VariationData& vd, int maxReadLength) {
    // Java filterAllSVStructures order: INV clusters, then DEL, then DUP. Each list is independent, so
    // order does not affect results; the INV lists must be collapsed so findINV sees mstart/mend/mlen.
    filterSVList(vd.svfinv3, maxReadLength);
    filterSVList(vd.svrinv3, maxReadLength);
    filterSVList(vd.svfinv5, maxReadLength);
    filterSVList(vd.svrinv5, maxReadLength);
    filterSVList(vd.svfdel, maxReadLength);
    filterSVList(vd.svrdel, maxReadLength);
    filterSVList(vd.svfdup, maxReadLength);
    filterSVList(vd.svrdup, maxReadLength);
}

void findDELdisc(VariationData& vd, Reference& ref, const Config& cfg, const Region& region, int maxReadLength) {
    (void)region;
    auto& NIV = vd.nonInsertionVariants;
    const int MINDIST = 8 * maxReadLength;

    auto buildTv = [](const Sclip& del) {
        Variation tv;
        tv.varsCount = 2 * del.varsCount; tv.highQualityReadsCount = 2 * del.varsCount;
        tv.varsCountOnForward = del.varsCount; tv.varsCountOnReverse = del.varsCount;
        tv.meanQuality = 2 * del.meanQuality; tv.meanPosition = 2 * del.meanPosition;
        tv.meanMappingQuality = 2 * del.meanMappingQuality; tv.numberOfMismatches = 2 * del.numberOfMismatches;
        return tv;
    };

    // forward clusters (svfdel): (svfdel) --> | ........ | <-- (svrdel)
    for (auto& del : vd.svfdel) {
        if (del.used) continue;
        if (del.varsCount < cfg.minReads + 5) continue;
        if (del.mstart <= del.end + MINDIST) continue;
        if (del.varsCount == 0 || del.meanMappingQuality / del.varsCount <= Config::DISCPAIRQUAL) continue;
        int mlen = del.mstart - del.end - maxReadLength / (del.varsCount + 1);
        if (!(mlen > 0 && mlen > MINDIST)) continue;
        int bp = del.end + (maxReadLength / (del.varsCount + 1)) / 2;
        if (del.softp != 0) bp = del.softp;
        ref.ensure(bp - 150, bp + 150);
        Variation& vref = getVariation(NIV, bp, "-" + std::to_string(mlen));
        vref.varsCount = 0;
        SVInfo& sv = vd.svInfoAt[bp];
        sv.type = "DEL";
        { auto it = vd.softClips3End.find(del.end + 1); if (it != vd.softClips3End.end()) sv.splits += it->second.varsCount; }
        { auto it = vd.softClips5End.find(del.mstart);  if (it != vd.softClips5End.end())  sv.splits += it->second.varsCount; }
        sv.pairs += del.varsCount;
        sv.clusters++;
        Variation tv = buildTv(del);
        adjCnt(vref, tv);
        if (!vd.refCoverage.count(bp)) vd.refCoverage[bp] = 2 * del.varsCount;
        del.used = true;
        markSVDel(del.end, del.mstart, vd.svrdel, maxReadLength);
    }

    // reverse clusters (svrdel)
    for (auto& del : vd.svrdel) {
        if (del.used) continue;
        if (del.varsCount < cfg.minReads + 5) continue;
        if (del.start <= del.mend + MINDIST) continue;
        if (del.varsCount == 0 || del.meanMappingQuality / del.varsCount <= Config::DISCPAIRQUAL) continue;
        int mlen = del.start - del.mend - maxReadLength / (del.varsCount + 1);
        if (!(mlen > 0 && mlen > MINDIST)) continue;
        int bp = del.mend + (maxReadLength / (del.varsCount + 1)) / 2;
        ref.ensure(bp - 150, bp + 150);
        Variation& vr = getVariation(NIV, bp, "-" + std::to_string(mlen));
        vr.varsCount = 0;
        SVInfo& sv = vd.svInfoAt[bp];
        sv.type = "DEL";
        { auto it = vd.softClips3End.find(del.mend + 1); if (it != vd.softClips3End.end()) sv.splits += it->second.varsCount; }
        { auto it = vd.softClips5End.find(del.start);   if (it != vd.softClips5End.end())  sv.splits += it->second.varsCount; }
        sv.pairs += del.varsCount;
        sv.clusters += 1;
        if (del.softp != 0) { auto it = vd.softClips5End.find(del.softp); if (it != vd.softClips5End.end()) it->second.used = true; }
        Variation tv = buildTv(del);
        adjCnt(vr, tv);
        if (!vd.refCoverage.count(bp)) vd.refCoverage[bp] = 2 * del.varsCount;
        if (vd.refCoverage.count(del.start) && vd.refCoverage[bp] < vd.refCoverage[del.start])
            vd.refCoverage[bp] = vd.refCoverage[del.start];
        del.used = true;
        ref.ensure(del.mstart - 100, del.mend + 100);
        markSVDel(del.mend, del.start, vd.svfdel, maxReadLength);
    }
}

// StructuralVariantsProcessor.findINVsub: pair-assisted inversion caller. Walks a discordant
// same-orientation INV cluster list (svfinv5/svrinv5/svfinv3/svrinv3), takes the cluster's dominant
// soft-clip position, finds the reciprocal breakpoint by matching the soft-clip consensus against the
// reverse-complemented reference (findMatchRev), and emits a "-<len>^<ins5><invM><ins3>" variation at
// the anchor position. The soft-clip reads are folded into the inversion and (dir==-1) subtracted from
// the reference variation via the 3-arg adjCnt, so the reference coverage is not double-counted.
// Returns after the first cluster that yields an inversion (mirrors the Java `return vref`).
static void findINVsub(std::vector<Sclip>& svref, int dir, int side,
                       VariationData& vd, Reference& ref, const Config& cfg,
                       const Region& region, int maxReadLength, const SVReloadFn& reload,
                       const std::vector<BamReader*>& bams) {
    auto& NIV = vd.nonInsertionVariants;
    for (Sclip& inv : svref) {
        if (inv.used) continue;
        if (inv.varsCount < cfg.minReads) continue;
        // dominant soft-clip position (max count; lowest position on tie), matching filterSV ordering.
        int softp = 0; { int bestc = -1; for (auto& [p, cnt] : inv.soft) if (cnt > bestc) { bestc = cnt; softp = p; } }
        std::map<int, Sclip>& sclip = dir == 1 ? vd.softClips3End : vd.softClips5End;

        // Pull the reciprocal breakpoint's reference AND reads into the window if it lies outside it
        // (Java's getReference(500) + partialPipeline reload over [mstart-200, mend+200]). The re-read
        // populates refCoverage at the far anchor, so a low-VAF inversion sitting on a high-coverage
        // locus is AF-filtered exactly as in VarDict (without it, its Depth collapses to the alt count
        // and it is emitted as a false positive).
        if (!(ref.has(inv.mstart) && ref.has(inv.mend))) {
            ref.ensure(inv.mstart - 500, inv.mend + 500);
            reload(inv.mstart, inv.mend);
        }

        int bp = 0; Sclip* scv = nullptr; std::string seq, extra;
        if (softp != 0) {
            auto it = sclip.find(softp); if (it == sclip.end()) continue;
            scv = &it->second; if (scv->used) continue;
            seq = findconseq(*scv); if (seq.empty()) continue;
            Match m = findMatchRev(seq, ref, softp, dir); bp = m.bp; extra = m.extra;
            if (bp == 0) { m = findMatchRev(seq, ref, softp, dir, Reference::SEED_2, 0); bp = m.bp; extra = m.extra; }
            if (bp == 0) continue;
        } else {
            int sp = dir == 1 ? inv.end : inv.start;
            for (int i = 1; i <= 2 * maxReadLength; i++) {
                int cp = sp + i * dir;
                auto it = sclip.find(cp); if (it == sclip.end()) continue;
                scv = &it->second; if (scv->used) continue;
                seq = findconseq(*scv); if (seq.empty()) continue;
                Match m = findMatchRev(seq, ref, cp, dir); bp = m.bp; extra = m.extra;
                if (bp == 0) { m = findMatchRev(seq, ref, cp, dir, Reference::SEED_2, 0); bp = m.bp; extra = m.extra; }
                if (bp == 0) continue;
                softp = cp;
                if ((dir == 1 && std::abs(bp - inv.mend) < Config::MINSVCDIST * maxReadLength)
                    || (dir == -1 && std::abs(bp - inv.mstart) < Config::MINSVCDIST * maxReadLength)) break;
            }
            if (bp == 0) continue;
        }
        if (side == 5) { if (dir == -1) bp--; }
        else { if (dir == 1) { bp++; if (bp != 0) softp--; } else { softp--; } }
        if (side == 3) { int tmp = bp; bp = softp; softp = tmp; }
        auto comp = [](char ch) { return complementBase(ch); };
        if ((dir == -1 && side == 5) || (dir == 1 && side == 3)) {
            while (ref.has(softp) && ref.has(bp) && ref.at(softp) == comp(ref.at(bp))) { softp++; if (softp != 0) bp--; }
        }
        while (ref.has(softp - 1) && ref.has(bp + 1) && ref.at(softp - 1) == comp(ref.at(bp + 1))) { softp--; if (softp != 0) bp++; }

        if (bp > softp && bp - softp > 150 && (bp - softp) / (double)std::abs(inv.mlen) < 1.5) {
            int len = bp - softp + 1;
            std::string ins5 = reverseComplement(joinRef(ref, bp - Config::SVFLANK + 1, bp));
            std::string ins3 = reverseComplement(joinRef(ref, softp, softp + Config::SVFLANK - 1));
            std::string ins = ins5 + "<inv" + std::to_string(len - 2 * Config::SVFLANK) + ">" + ins3;
            if (len - 2 * Config::SVFLANK <= 0) ins = reverseComplement(joinRef(ref, softp, bp));
            if (dir == 1 && !extra.empty()) { extra = reverseComplement(extra); ins = extra + ins; }
            else if (dir == -1 && !extra.empty()) { ins = ins + extra; }
            std::string gt = "-" + std::to_string(len) + "^" + ins;

            Variation& vref = getVariation(NIV, softp, gt);
            inv.used = true; vref.pstd = true; vref.qstd = true;
            SVInfo& sv = vd.svInfoAt[softp];
            sv.type = "INV"; sv.splits += scv->varsCount; sv.pairs += inv.varsCount; sv.clusters++;

            Variation* vrefSoftp = (dir == -1 && ref.has(softp))
                ? getVariationMaybe(NIV, softp, ref.at(softp)) : nullptr;
            adjCnt(vref, *scv, vrefSoftp);
            vd.refCoverage[softp] = vd.refCoverage.count(softp - 1) ? vd.refCoverage[softp - 1] : inv.varsCount;
            scv->used = true;
            // Java re-runs realigndel on the single {softp:{gt:inv.varsCount}} INV "deletion" hash
            // (StructuralVariantsProcessor.findINVsub l.658-667). The opposite-strand soft clip at softp
            // (reverse-complement of scv) gets attracted into the inversion, raising AltDepth/coverage.
            realignOneDel(vd, ref, cfg, region, maxReadLength, bams, softp, gt, inv.varsCount);
            return;
        }
    }
}

void findINV(VariationData& vd, Reference& ref, const Config& cfg, const Region& region,
             int maxReadLength, const SVReloadFn& reload, const std::vector<BamReader*>& bams) {
    if (cfg.disableSV) return;
    findINVsub(vd.svfinv5, 1, 5, vd, ref, cfg, region, maxReadLength, reload, bams);
    findINVsub(vd.svrinv5, -1, 5, vd, ref, cfg, region, maxReadLength, reload, bams);
    findINVsub(vd.svfinv3, 1, 3, vd, ref, cfg, region, maxReadLength, reload, bams);
    findINVsub(vd.svrinv3, -1, 3, vd, ref, cfg, region, maxReadLength, reload, bams);
}

// StructuralVariantsProcessor.findsv: split-read SVs on 3'/5' soft clips. The forward findMatch branch
// (candidate DEL/DUP) needs discordant-pair support via checkPairs, which is unported (no discordant SV
// pairs are collected on WES targets, so checkPairs would return 0 and Java just `continue`s); the DUP
// sub-branch is empty in Java. Only the candidate-inversion path (findMatch fails, findMatchRev matches
// the reverse strand) emits variations here, producing the split-read <INV> calls (SVinfo cnt-0-0).
void findsv(VariationData& vd, Reference& ref, const Config& cfg, const Region& region, int maxReadLength) {
    (void)maxReadLength;
    auto& NIV = vd.nonInsertionVariants;

    // fillAndSortTmpSV: unused soft clips within the current segment, sorted by count descending.
    struct SortSclip { int position; Sclip* sc; int count; };
    auto fillAndSort = [&](std::map<int, Sclip>& clips) {
        std::vector<SortSclip> tmp;
        for (auto& [pos, sc] : clips) {
            if (sc.used) continue;
            if (pos < region.start || pos > region.end) continue;
            tmp.push_back({pos, &sc, sc.varsCount});
        }
        std::stable_sort(tmp.begin(), tmp.end(),
                         [](const SortSclip& a, const SortSclip& b) { return a.count > b.count; });
        return tmp;
    };
    auto& refseq = ref;   // reference.referenceSequences accessor
    auto isHasAndEq = [&](char c, int pos) { return refseq.has(pos) && refseq.at(pos) == c; };

    // 5' soft clips
    for (auto& t5 : fillAndSort(vd.softClips5End)) {
        int p5 = t5.position;
        Sclip& sc5v = *t5.sc;
        int cnt5 = t5.count;
        if (cnt5 < cfg.minReads) break;
        if (sc5v.used) continue;
        std::string seq = findconseq(sc5v);
        if (seq.empty() || (int)seq.size() < Reference::SEED_2) continue;
        Match match = findMatch(seq, ref, p5, -1, Reference::SEED_1, 3);
        int bp = match.bp;
        std::string EXTRA = match.extra;
        if (bp != 0) continue;   // candidate DEL/DUP (checkPairs unported -> no output; DUP is empty)
        // candidate inversion
        Match matchRev = findMatchRev(seq, ref, p5, -1);
        bp = matchRev.bp; EXTRA = matchRev.extra;
        if (bp == 0) continue;
        if (!(std::abs(bp - p5) > Config::SVFLANK)) continue;
        if (bp <= p5) { int temp = bp; bp = p5; p5 = temp; }
        bp--;
        while (refseq.has(bp + 1) && isHasAndEq(complementBase(refseq.at(bp + 1)), p5 - 1)) {
            p5--; if (p5 != 0) bp++;
        }
        std::string ins5 = reverseComplement(joinRef(ref, bp - Config::SVFLANK + 1, bp));
        std::string ins3 = reverseComplement(joinRef(ref, p5, p5 + Config::SVFLANK - 1));
        int mid = bp - p5 - (int)ins5.size() - (int)ins3.size() + 1;
        std::string vn = "-" + std::to_string(bp - p5 + 1) + "^" + ins5 + "<inv" + std::to_string(mid) + ">" + ins3 + EXTRA;
        if (mid <= 0) {
            std::string tins = reverseComplement(joinRef(ref, p5, bp));
            vn = "-" + std::to_string(bp - p5 + 1) + "^" + tins + EXTRA;
        }
        Variation& vref = getVariation(NIV, p5, vn);
        SVInfo& sv = vd.svInfoAt[p5];
        sv.type = "INV"; sv.splits += cnt5;
        adjCnt(vref, sc5v);
        vd.refCoverage[p5] += cnt5;
        if (vd.refCoverage.count(bp) && vd.refCoverage[p5] < vd.refCoverage[bp]) vd.refCoverage[p5] = vd.refCoverage[bp];
    }

    // 3' soft clips
    for (auto& t3 : fillAndSort(vd.softClips3End)) {
        int p3 = t3.position;
        Sclip& sc3v = *t3.sc;
        int cnt3 = t3.count;
        if (cnt3 < cfg.minReads) break;
        if (sc3v.used) continue;
        std::string seq = findconseq(sc3v);
        if (seq.empty() || (int)seq.size() < Reference::SEED_2) continue;
        Match match = findMatch(seq, ref, p3, 1, Reference::SEED_1, 3);
        int bp = match.bp;
        std::string EXTRA = match.extra;
        if (bp != 0) continue;   // candidate DEL/DUP (checkPairs unported -> no output; DUP is empty)
        // candidate inversion
        Match matchRev = findMatchRev(seq, ref, p3, 1);
        bp = matchRev.bp; EXTRA = matchRev.extra;
        if (bp == 0) continue;
        if (std::abs(bp - p3) <= Config::SVFLANK) continue;
        if (bp < p3) { int tmp = bp; bp = p3; p3 = tmp; p3++; bp--; }
        while (refseq.has(bp + 1) && isHasAndEq(complementBase(refseq.at(bp + 1)), p3 - 1)) {
            p3--; if (p3 != 0) bp++;
        }
        std::string ins5 = reverseComplement(joinRef(ref, bp - Config::SVFLANK + 1, bp));
        std::string ins3 = reverseComplement(joinRef(ref, p3, p3 + Config::SVFLANK - 1));
        int mid = bp - p3 - 2 * Config::SVFLANK + 1;
        std::string vn = "-" + std::to_string(bp - p3 + 1) + "^" + EXTRA + ins5 + "<inv" + std::to_string(mid) + ">" + ins3;
        if (mid <= 0) {
            std::string tins = reverseComplement(joinRef(ref, p3, bp));
            vn = "-" + std::to_string(bp - p3 + 1) + "^" + EXTRA + tins;
        }
        Variation& vref = getVariation(NIV, p3, vn);
        SVInfo& sv = vd.svInfoAt[p3];
        sv.type = "INV"; sv.splits += cnt3;
        adjCnt(vref, sc3v);
        vd.refCoverage[p3] += cnt3;
        if (vd.refCoverage.count(bp) && vd.refCoverage[p3] < vd.refCoverage[bp]) vd.refCoverage[p3] = vd.refCoverage[bp];
    }
}

} // namespace vardict
