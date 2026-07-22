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

static Variation& getVariation(std::map<int, VarMap>& hash, int pos, const std::string& key) {
    return hash[pos][key];
}
static Variation* getVariationMaybe(std::map<int, VarMap>& hash, int pos, char refBase) {
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
    if (!SEQ.empty() && islowcomplexseq(SEQ)) sc.used = true;
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

void adjustMNP(VariationData& vd, Reference&, const Config&, const Region&) {
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
    }
}

// ---- realignins ---------------------------------------------------------------------------------

void realignins(VariationData& vd, Reference& ref, const Config& cfg, const Region& region, int maxReadLength) {
    struct Item { int position; std::string desc; int count; };
    std::vector<Item> tmp;
    for (auto& [pos, m] : vd.positionToInsertionCount) for (auto& [d, c] : m) tmp.push_back({pos, d, c});
    std::sort(tmp.begin(), tmp.end(), [](const Item& a, const Item& b) {
        if (a.count != b.count) return a.count > b.count;
        if (a.position != b.position) return a.position < b.position;
        return a.desc > b.desc;
    });
    auto& NIV = vd.nonInsertionVariants;
    for (const auto& t : tmp) {
        int position = t.position; const std::string& vn = t.desc; int insertionCount = t.count;
        // BEGIN_PLUS_ATGC: insert = sequence after '+'. (dup/&/#/^ grammar not produced by this port.)
        if (vn.empty() || vn[0] != '+') continue;
        std::string insert = vn.substr(1);
        for (auto& c : insert) if (!(c=='A'||c=='C'||c=='G'||c=='T')) { insert.clear(); break; }
        if (insert.empty()) continue;
        std::string extra, compm; int newdel = 0; std::string ins3;
        int inslen = (int)insert.size();
        std::string tn = insert;   // vn with +,&,#,^N,^ stripped == the plain inserted seq

        int wustart = position - 150 > 1 ? position - 150 : 1;
        std::string wupseq = joinRef(ref, wustart, position) + tn;
        int sanend = position + (int)vn.size() + 100;
        std::string sanpseq = tn + joinRef(ref, position + (int)extra.size() + 1 + (int)compm.size() + newdel, sanend);
        MismatchResult f3 = findMM3(vd, ref, position + 1, sanpseq);
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

void realigndel(VariationData& vd, Reference& ref, const Config& cfg, const Region& region, int maxReadLength) {
    struct Item { int position; std::string desc; int count; };
    std::vector<Item> tmp;
    for (auto& [pos, m] : vd.positionToDeletionCount) for (auto& [d, c] : m) tmp.push_back({pos, d, c});
    std::sort(tmp.begin(), tmp.end(), [](const Item& a, const Item& b) {
        if (a.count != b.count) return a.count > b.count;
        if (a.position != b.position) return a.position < b.position;
        return a.desc > b.desc;
    });
    auto& NIV = vd.nonInsertionVariants;
    for (const auto& t : tmp) {
        int p = t.position; const std::string& vn = t.desc; int dcnt = t.count;
        // BEGIN_MINUS_NUMBER: dellen. (^N / SV / & grammar not produced by this port.)
        if (vn.empty() || vn[0] != '-') continue;
        int dellen = std::atoi(vn.c_str() + 1);
        std::string extra, extrains;
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
    }
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

// ---- realignlgdel (large deletions from soft-clip breakpoints) -----------------------------------
// Faithful port of the findbp path of VariationRealigner.realignlgdel. The bp==0 fallback (seed-based
// findMatch + discordant-pair SV clusters + partialPipeline on extended regions) is gated off — those
// require the SV subsystem not yet ported. SV output markers are likewise skipped (no SV column yet).

void realignlgdel(VariationData& vd, Reference& ref, const Config& cfg, const Region& region, int maxReadLength) {
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
        if (bp == 0) continue; // seed/SV fallback not ported
        int dellen = p - bp;
        std::string extra;
        std::string gt = "-" + std::to_string(dellen);
        int en = 0;
        while (en < (int)seq.size() && !(ref.has(bp - en - 1) && seq[en] == ref.at(bp - en - 1))) {
            extra += seq[en]; en++;
        }
        if (!extra.empty()) { extra = reverseStr(extra); gt = "-" + std::to_string(dellen) + "&" + extra; bp -= (int)extra.size(); }
        // 3' breakpoint (sc3p)
        int n = 0;
        if (extra.empty()) {
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
        if (bp == 0) continue;
        int dellen = bp - p;
        std::string extra; int en = 0;
        while (en < (int)seq.size() && !(ref.has(bp + en) && seq[en] == ref.at(bp + en))) { extra += seq[en]; en++; }
        std::string gt = "-" + std::to_string(dellen);
        bp = p; // set to 5'
        if (!extra.empty()) gt = "-" + std::to_string(dellen) + "&" + extra;
        else {
            while (ref.has(bp - 1) && ref.has(bp + dellen - 1) && ref.at(bp - 1) == ref.at(bp + dellen - 1)) bp--;
        }
        Variation& tv = getVariation(NIV, bp, gt);
        tv.qstd = true; tv.pstd = true;
        if (dellen < cfg.indelsize)
            for (int tp = bp; tp < bp + dellen + (int)extra.size(); ++tp) vd.refCoverage[tp] += sc3v.varsCount;
        if (!vd.refCoverage.count(bp)) vd.refCoverage[bp] = vd.refCoverage.count(p - 1) ? vd.refCoverage[p - 1] : sc3v.varsCount;
        sc3v.meanPosition += (double)dellen * sc3v.varsCount;
        adjCnt(tv, sc3v); sc3v.used = true;
    }
}

} // namespace vardict
