#include "cigar_modifier.hpp"
#include "util.hpp"

namespace vardict {

static inline bool refEq(Reference& r, int rp, const std::string& s, int sp) {
    return r.has(rp) && sp >= 0 && sp < (int)s.size() && r.at(rp) == s[sp];
}
static inline bool refNeq(Reference& r, int rp, const std::string& s, int sp) {
    return r.has(rp) && sp >= 0 && sp < (int)s.size() && r.at(rp) != s[sp];
}
static int countIndel(const Cig& c) { int n = 0; for (auto& e : c) if (e.second == 'I' || e.second == 'D') n++; return n; }
static bool isM(char op) { return op == 'M' || op == '=' || op == 'X'; }
static bool consumesRef(char op) { return isM(op) || op == 'D' || op == 'N'; }
static bool consumesRead(char op) { return isM(op) || op == 'I' || op == 'S'; }

// captureMisSoftlyMS: trailing "..M##S" -> extend the match into the trailing soft-clip where bases match.
static void captureMisSoftlyMS(int position, Cig& cig, const std::string& seq, const std::vector<int>& qual, Reference& ref) {
    int k = (int)cig.size();
    if (k < 2 || cig[k - 1].second != 'S' || !isM(cig[k - 2].second)) return;
    int mch = cig[k - 2].first, soft = cig[k - 1].first;
    int refoff = position, rdoff = 0;
    for (int i = 0; i < k - 1; ++i) { if (consumesRef(cig[i].second)) refoff += cig[i].first; if (consumesRead(cig[i].second)) rdoff += cig[i].first; }
    int rn = 0;
    while (rn < soft && refEq(ref, refoff + rn, seq, rdoff + rn) && qual[rdoff + rn] > Config::LOWQUAL) rn++;
    if (rn > 0) { mch += rn; soft -= rn; }
    if (soft > 0) {
        int rn2 = 0; std::string RN;
        while (rn2 + 1 < soft && refEq(ref, refoff + rn + rn2 + 1, seq, rdoff + rn + rn2 + 1) && qual[rdoff + rn + rn2 + 1] > Config::LOWQUAL) {
            rn2++; if (ref.has(refoff + rn + rn2 + 1)) { char c = ref.at(refoff + rn + rn2 + 1); if (RN.find(c) == std::string::npos) RN += c; }
        }
        if (rn2 > 4 && (int)RN.size() > 1) { mch += rn2 + 1; soft -= rn2 + 1; }
    }
    cig[k - 2].first = mch;
    if (soft > 0) cig[k - 1].first = soft; else cig.pop_back();
}

// captureMisSoftly3Mismatches: trailing "..##M" (no trailing S) -> soft-clip up to 3 trailing mismatched bases.
static void captureMisSoftly3Mismatches(int position, Cig& cig, const std::string& seq, Reference& ref) {
    int k = (int)cig.size();
    if (k < 1 || !isM(cig[k - 1].second)) return;
    int mch = cig[k - 1].first;
    int refoff = position, rdoff = 0;
    for (auto& e : cig) { if (consumesRef(e.second)) refoff += e.first; if (consumesRead(e.second)) rdoff += e.first; }
    int rn = 0, rrn = 0, rmch = 0;
    while (rrn < mch && rn < mch) {
        if (!ref.has(refoff - rrn - 1)) break;
        if (rrn < rdoff && refNeq(ref, refoff - rrn - 1, seq, rdoff - rrn - 1)) { rn = rrn + 1; rmch = 0; }
        else if (rrn < rdoff && refEq(ref, refoff - rrn - 1, seq, rdoff - rrn - 1)) rmch++;
        rrn++;
        if (rmch >= 3) break;
    }
    if (rn > 0 && rn <= 3) { cig[k - 1].first = mch - rn; cig.push_back({rn, 'S'}); }
}

// combineDigSDigM: leading "##S##M" -> extend match into the leading soft-clip / soft-clip leading mismatches.
static void combineDigSDigM(int& position, Cig& cig, const std::string& seq, const std::vector<int>& qual, Reference& ref) {
    if (cig.size() < 2 || cig[0].second != 'S' || !isM(cig[1].second)) return;
    int soft = cig[0].first, mch = cig[1].first;
    int rn = 0;
    while (rn < soft && refEq(ref, position - rn - 1, seq, soft - rn - 1) && qual[soft - rn - 1] > Config::LOWQUAL) rn++;
    if (rn > 0) { mch += rn; soft -= rn; position -= rn; rn = 0; }
    if (soft > 0) {
        int rn3 = 0; std::string RN;
        while (rn3 + 1 < soft && refEq(ref, position - rn3 - 2, seq, soft - rn3 - 2) && qual[soft - rn3 - 2] > Config::LOWQUAL) {
            rn3++; if (ref.has(position - rn3 - 2)) { char c = ref.at(position - rn3 - 2); if (RN.find(c) == std::string::npos) RN += c; }
        }
        if ((rn3 > 4 && (int)RN.size() > 1) || refEq(ref, position - 1, seq, soft - 1)) { mch += rn3 + 1; soft -= rn3 + 1; position -= rn3 + 1; }
        if (rn3 == 0) {
            int rrn = 0, rmch = 0, r2 = 0;
            while (rrn < mch && r2 < mch) {
                if (!ref.has(position + rrn)) break;
                if (refNeq(ref, position + rrn, seq, soft + rrn)) { r2 = rrn + 1; rmch = 0; }
                else if (refEq(ref, position + rrn, seq, soft + rrn)) rmch++;
                rrn++;
                if (rmch >= 3) break;
            }
            if (r2 > 0 && r2 < mch) { soft += r2; mch -= r2; position += r2; }
        }
    }
    if (soft > 0) { cig[0].first = soft; cig[1].first = mch; }
    else { cig.erase(cig.begin()); cig[0] = {mch, 'M'}; }
}

// combineBeginDigM: leading "##M" (no leading S) -> soft-clip up to 3 leading mismatched bases.
static void combineBeginDigM(int& position, Cig& cig, const std::string& seq, Reference& ref) {
    if (cig.empty() || !isM(cig[0].second)) return;
    int mch = cig[0].first, rn = 0, rrn = 0, rmch = 0;
    while (rrn < mch && rn < mch) {
        if (!ref.has(position + rrn)) break;
        if (refNeq(ref, position + rrn, seq, rrn)) { rn = rrn + 1; rmch = 0; }
        else if (refEq(ref, position + rrn, seq, rrn)) rmch++;
        rrn++;
        if (rmch >= 3) break;
    }
    if (rn > 0 && rn <= 3) { cig[0].first = mch - rn; cig.insert(cig.begin(), {rn, 'S'}); position += rn; }
}

void modifyCigar(int& position, Cig& cig, std::string& seq, std::vector<int>& qual,
                 Reference& ref, int maxReadLength, const Config& cfg) {
    if (cig.empty()) return;
    if (cig.front().second == 'D') { position += cig.front().first; cig.erase(cig.begin()); }
    if (!cig.empty() && cig.back().second == 'D') cig.pop_back();
    if (cig.empty()) return;
    if (cig.front().second == 'I') cig.front().second = 'S';
    if (cig.back().second == 'I') cig.back().second = 'S';

    // chimeric-seed clip removal (SEED_2): drop a large leading/trailing soft-clip whose reverse-
    // complement maps uniquely near the read.
    if (!cfg.chimeric) {
        if (cig.front().second == 'S' && cig.front().first >= Reference::SEED_2) {
            int el = cig.front().first;
            const auto* pos = ref.seedPositions(reverseComplement(seq.substr(0, el)).substr(0, Reference::SEED_2));
            if (pos && pos->size() == 1 && std::abs(position - (*pos)[0]) < 2 * maxReadLength) {
                cig.erase(cig.begin()); seq = seq.substr(el); qual.erase(qual.begin(), qual.begin() + el);
            }
        } else if (cig.back().second == 'S' && cig.back().first >= Reference::SEED_2) {
            int el = cig.back().first;
            std::string rc = reverseComplement(seq.substr(seq.size() - el, el));
            const auto* pos = ref.seedPositions(rc.substr(rc.size() - Reference::SEED_2, Reference::SEED_2));
            if (pos && pos->size() == 1 && std::abs(position - (*pos)[0]) < 2 * maxReadLength) {
                cig.pop_back(); seq = seq.substr(0, seq.size() - el); qual.resize(qual.size() - el);
            }
        }
    }

    // CigarModifier leading indel-normalization (partial port of the while(flag && indel>0) loop):
    // ^(\d+)S(\d+)M(\d+)([ID]) with the anchored match <= 10 bp. The short match wedged between a
    // soft-clip and an indel is an unreliable anchor, so VarDict folds soft-clip + match (+ any
    // inserted bases) into a single soft-clip and advances the reference position past the consumed
    // match (+ any deleted bases); its realigner re-finds the indel from the clip if it is real.
    // (BEGIN_NUMBER_S_NUMBER_M_NUMBER_IorD, CigarModifier.java:156)
    if (cig.size() >= 3 && cig[0].second == 'S' && cig[1].second == 'M' && cig[1].first <= 10 &&
        (cig[2].second == 'I' || cig[2].second == 'D')) {
        int s = cig[0].first, m = cig[1].first, d = cig[2].first;
        bool isI = (cig[2].second == 'I');
        position += m + (isI ? 0 : d);
        cig.erase(cig.begin(), cig.begin() + 3);
        cig.insert(cig.begin(), {s + m + (isI ? d : 0), 'S'});
    }

    // Reads with (remaining) indels are left to the (un-ported) indel-collapse loop.
    if (countIndel(cig) > 0) return;

    // Trailing: ..M##S -> captureMisSoftlyMS ; else ..##M -> captureMisSoftly3Mismatches
    if (cig.size() >= 2 && cig.back().second == 'S' && isM(cig[cig.size() - 2].second))
        captureMisSoftlyMS(position, cig, seq, qual, ref);
    else if (!cig.empty() && isM(cig.back().second))
        captureMisSoftly3Mismatches(position, cig, seq, ref);

    // Leading: ##S##M -> combineDigSDigM ; else ##M -> combineBeginDigM
    if (cig.size() >= 2 && cig[0].second == 'S' && isM(cig[1].second))
        combineDigSDigM(position, cig, seq, qual, ref);
    else if (!cig.empty() && isM(cig[0].second))
        combineBeginDigM(position, cig, seq, ref);
}

} // namespace vardict
