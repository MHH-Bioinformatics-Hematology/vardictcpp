#include "cigar_modifier.hpp"
#include "util.hpp"
#include <regex>

namespace vardict {

static inline bool refEq(Reference& r, int rp, const std::string& s, int sp) {
    return r.has(rp) && sp >= 0 && sp < (int)s.size() && r.at(rp) == s[sp];
}
static inline bool refNeq(Reference& r, int rp, const std::string& s, int sp) {
    return r.has(rp) && sp >= 0 && sp < (int)s.size() && r.at(rp) != s[sp];
}
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
    // Java resets rn to 0 here (captureMisSoftlyMS line 430) and the second scan REUSES rn (not a
    // fresh counter): it therefore rescans from refoff+1, re-testing the base right after the M/S
    // boundary. Carrying the first scan's rn into the index skips past that base and over-extends.
    if (rn > 0) { mch += rn; soft -= rn; rn = 0; }
    if (soft > 0) {
        std::string RN;
        while (rn + 1 < soft && refEq(ref, refoff + rn + 1, seq, rdoff + rn + 1) && qual[rdoff + rn + 1] > Config::LOWQUAL) {
            rn++; if (ref.has(refoff + rn + 1)) { char c = ref.at(refoff + rn + 1); if (RN.find(c) == std::string::npos) RN += c; }
        }
        if (rn > 4 && (int)RN.size() > 1) { mch += rn + 1; soft -= rn + 1; }
        // Java captureMisSoftlyMS `if (rn == 0)` block: when the forward scan found no match,
        // walk backward from the M/S boundary and soft-clip back to the last mismatch found
        // within a run of <3 consecutive matches (moves the soft-clip start earlier).
        if (rn == 0) {
            int rrn = 0, rmch = 0, rnb = 0;
            while (rrn < mch && rnb < mch) {
                if (!ref.has(refoff - rrn - 1)) break;
                if (rrn < rdoff && refNeq(ref, refoff - rrn - 1, seq, rdoff - rrn - 1)) { rnb = rrn + 1; rmch = 0; }
                else if (rrn < rdoff && refEq(ref, refoff - rrn - 1, seq, rdoff - rrn - 1)) rmch++;
                rrn++;
                if (rmch >= 3) break;
            }
            if (rnb > 0 && rnb < mch) { soft += rnb; mch -= rnb; }
        }
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

// ---- Faithful string-CIGAR port of CigarModifier's `while (flag && indel > 0)` indel-collapse
//      loop (CigarModifier.java:132-227). Only mutates the CIGAR string and `position`; the read
//      sequence/quality are never touched by these transforms, so we operate on a string and
//      convert back to the Cig vector afterwards. Runs only for reads that carry an indel, so the
//      per-read cost (short CIGAR, a handful of regex probes) stays off the common no-indel path.

static std::string cigToStr(const Cig& c) {
    std::string s;
    for (auto& e : c) { s += std::to_string(e.first); s += e.second; }
    return s;
}
static Cig strToCig(const std::string& s) {
    Cig c; int i = 0, n = (int)s.size();
    while (i < n) {
        int num = 0;
        while (i < n && s[i] >= '0' && s[i] <= '9') { num = num * 10 + (s[i] - '0'); i++; }
        if (i < n) { c.push_back({num, s[i]}); i++; }
    }
    return c;
}
// Utils.sum(globalFind(pattern, s)): sum of group(1) integers over all non-overlapping matches.
static int sumGroup1(const std::string& s, const std::regex& re) {
    int total = 0;
    for (auto it = std::sregex_iterator(s.begin(), s.end(), re), end = std::sregex_iterator(); it != end; ++it)
        total += std::stoi((*it)[1].str());
    return total;
}

static void indelCollapseLoop(std::string& cs, int& position, const std::string& seq, Reference& ref) {
    static const std::regex BEGIN_S_ID(R"(^(\d+)S(\d+)([ID]))");
    static const std::regex END_ID_S(R"((\d+)([ID])(\d+)S$)");
    static const std::regex BEGIN_S_M_ID(R"(^(\d+)S(\d+)M(\d+)([ID]))");
    static const std::regex END_ID_M_S(R"((\d+)([ID])(\d+)M(\d+)S$)");
    static const std::regex BEGIN_DM_ID_M(R"(^(\d)M(\d+)([ID])(\d+)M)");
    static const std::regex END_ID_DM(R"((\d+)([ID])(\d)M$)");
    static const std::regex TWO_DEL_INS(R"(^(.*?)(\d+)M(\d+)D(\d+)M(\d+)I(\d+)M(\d+)D(\d+)M)");
    static const std::regex THREE_DEL(R"(^(.*?)(\d+)M(\d+)D(\d+)M(\d+)D(\d+)M(\d+)D(\d+)M)");
    static const std::regex THREE_INDEL(R"(^(.*?)(\d+)M(\d+)([DI])(\d+)M(\d+)([DI])(\d+)M(\d+)([DI])(\d+)M)");
    static const std::regex DIG_D_M_DI_I(R"((\d+)D(\d+)M(\d+)([DI])(\d+I)?)");
    static const std::regex DIG_I_M_DI_I(R"((\d+)I(\d+)M(\d+)([DI])(\d+I)?)");
    static const std::regex NOTDIG_I_M_DI_I(R"((\D)(\d+)I(\d+)M(\d+)([DI])(\d+I)?)");
    static const std::regex DIG_D_D(R"((\d+)D(\d+)D)");
    static const std::regex DIG_I_I(R"((\d+)I(\d+)I)");
    static const std::regex NUM_MND(R"((\d+)[MND])");
    static const std::regex NUM_MIS(R"((\d+)[MIS])");

    bool flag = true;
    std::smatch m;
    while (flag) {
        flag = false;
        // ^(\d+)S(\d+)([ID]): fold a leading soft-clip + indel into one soft-clip (D advances position).
        if (std::regex_search(cs, m, BEGIN_S_ID)) {
            int g1 = std::stoi(m[1]), g2 = std::stoi(m[2]); char op = m[3].str()[0];
            int tslen = g1 + (op == 'I' ? g2 : 0);
            if (op == 'D') position += g2;
            cs = std::to_string(tslen) + "S" + m.suffix().str();
            flag = true;
        }
        // (\d+)([ID])(\d+)S$: fold a trailing indel + soft-clip into one soft-clip.
        if (std::regex_search(cs, m, END_ID_S)) {
            int g1 = std::stoi(m[1]); char op = m[2].str()[0]; int g3 = std::stoi(m[3]);
            int tslen = g3 + (op == 'I' ? g1 : 0);
            cs = m.prefix().str() + std::to_string(tslen) + "S";
            flag = true;
        }
        // ^(\d+)S(\d+)M(\d+)([ID]) with M<=10: fold clip+short match+indel into one clip.
        if (std::regex_search(cs, m, BEGIN_S_M_ID)) {
            int g1 = std::stoi(m[1]), g2 = std::stoi(m[2]), g3 = std::stoi(m[3]); char op = m[4].str()[0];
            if (g2 <= 10) {
                int tslen = g1 + g2 + (op == 'I' ? g3 : 0);
                position += g2 + (op == 'D' ? g3 : 0);
                cs = std::to_string(tslen) + "S" + m.suffix().str();
                flag = true;
            }
        }
        // (\d+)([ID])(\d+)M(\d+)S$ with M<=10: fold indel+short match+clip into one clip.
        if (std::regex_search(cs, m, END_ID_M_S)) {
            int g1 = std::stoi(m[1]); char op = m[2].str()[0]; int g3 = std::stoi(m[3]), g4 = std::stoi(m[4]);
            if (g3 <= 10) {
                int tslen = g4 + g3 + (op == 'I' ? g1 : 0);
                cs = m.prefix().str() + std::to_string(tslen) + "S";
                flag = true;
            }
        }
        // ^(\d)M(\d+)([ID])(\d+)M: beginDigitMNumberIorDNumberM -- 1-9bp leading match before an indel;
        // fold into a soft-clip extended over leading mismatches of the following match.
        if (std::regex_search(cs, m, BEGIN_DM_ID_M)) {
            int tmid = std::stoi(m[1]), g2 = std::stoi(m[2]); char op = m[3].str()[0]; int mlen = std::stoi(m[4]);
            int tslen = tmid + (op == 'I' ? g2 : 0);
            position += tmid + (op == 'D' ? g2 : 0);
            int tn = 0;
            while (tn < mlen && refNeq(ref, position + tn, seq, tslen + tn)) tn++;
            tslen += tn; mlen -= tn; position += tn;
            cs = std::to_string(tslen) + "S" + std::to_string(mlen) + "M" + m.suffix().str();
            flag = true;
        }
        // (\d+)([ID])(\d)M$: trailing indel + 1-9bp match -> single trailing soft-clip.
        if (std::regex_search(cs, m, END_ID_DM)) {
            int g1 = std::stoi(m[1]); char op = m[2].str()[0]; int tmid = std::stoi(m[3]);
            int tslen = tmid + (op == 'I' ? g1 : 0);
            cs = m.prefix().str() + std::to_string(tslen) + "S";
            flag = true;
        }
        // Three close indels -> one D/I complex (else-if chain, exactly as Java).
        if (std::regex_search(cs, m, TWO_DEL_INS)) {              // twoDeletionsInsertionToComplex
            std::string ov5 = m[1].str();
            int g2 = std::stoi(m[2]), g3 = std::stoi(m[3]), g4 = std::stoi(m[4]), g5 = std::stoi(m[5]),
                g6 = std::stoi(m[6]), g7 = std::stoi(m[7]), g8 = std::stoi(m[8]);
            int tslen = g4 + g5 + g6, dlen = g3 + g4 + g6 + g7, mid = g4 + g6;
            int refoff = position + g2, rdoff = g2, RDOFF = g2, rm = g8;
            if (!ov5.empty()) { refoff += sumGroup1(ov5, NUM_MND); rdoff += sumGroup1(ov5, NUM_MIS); }
            int rn = 0; while (rdoff + rn < (int)seq.size() && refEq(ref, refoff + rn, seq, rdoff + rn)) rn++;
            RDOFF += rn; dlen -= rn; tslen -= rn;
            std::string nc = std::to_string(RDOFF) + "M";
            if (tslen <= 0) { dlen -= tslen; rm += tslen; nc += std::to_string(dlen) + "D" + std::to_string(rm) + "M"; }
            else { nc += std::to_string(dlen) + "D" + std::to_string(tslen) + "I" + std::to_string(rm) + "M"; }
            if (mid <= 15) { cs = ov5 + nc + m.suffix().str(); flag = true; }
        } else if (std::regex_search(cs, m, THREE_DEL)) {          // threeDeletions
            std::string ov5 = m[1].str();
            int g2 = std::stoi(m[2]), g3 = std::stoi(m[3]), g4 = std::stoi(m[4]), g5 = std::stoi(m[5]),
                g6 = std::stoi(m[6]), g7 = std::stoi(m[7]), g8 = std::stoi(m[8]);
            int tslen = g4 + g6, dlen = g3 + g4 + g5 + g6 + g7, mid = g4 + g6;
            int refoff = position + g2, rdoff = g2, RDOFF = g2, rm = g8;
            if (!ov5.empty()) { refoff += sumGroup1(ov5, NUM_MND); rdoff += sumGroup1(ov5, NUM_MIS); }
            int rn = 0; while (rdoff + rn < (int)seq.size() && refEq(ref, refoff + rn, seq, rdoff + rn)) rn++;
            RDOFF += rn; dlen -= rn; tslen -= rn;
            std::string nc = std::to_string(RDOFF) + "M";
            if (tslen <= 0) { dlen -= tslen; rm += tslen; nc += std::to_string(dlen) + "D" + std::to_string(rm) + "M"; }
            else { nc += std::to_string(dlen) + "D" + std::to_string(tslen) + "I" + std::to_string(rm) + "M"; }
            if (mid <= 15) { cs = ov5 + nc + m.suffix().str(); flag = true; }
        } else if (std::regex_search(cs, m, THREE_INDEL)) {        // threeIndels
            std::string ov5 = m[1].str();
            int g2 = std::stoi(m[2]), g3 = std::stoi(m[3]); char o1 = m[4].str()[0];
            int g5 = std::stoi(m[5]), g6 = std::stoi(m[6]); char o2 = m[7].str()[0];
            int g8 = std::stoi(m[8]), g9 = std::stoi(m[9]); char o3 = m[10].str()[0];
            int g11 = std::stoi(m[11]);
            int tslen = g5 + g8; if (o1 == 'I') tslen += g3; if (o2 == 'I') tslen += g6; if (o3 == 'I') tslen += g9;
            int dlen = g5 + g8;  if (o1 == 'D') dlen += g3;  if (o2 == 'D') dlen += g6;  if (o3 == 'D') dlen += g9;
            int mid = g5 + g8;
            int refoff = position + g2, rdoff = g2, RDOFF = g2, rm = g11;
            if (!ov5.empty()) { refoff += sumGroup1(ov5, NUM_MND); rdoff += sumGroup1(ov5, NUM_MIS); }
            int rn = 0; while (rdoff + rn < (int)seq.size() && refEq(ref, refoff + rn, seq, rdoff + rn)) rn++;
            RDOFF += rn; dlen -= rn; tslen -= rn;
            std::string nc = std::to_string(RDOFF) + "M";
            if (tslen <= 0) {
                dlen -= tslen; rm += tslen;
                if (dlen == 0) { RDOFF = RDOFF + rm; nc = std::to_string(RDOFF) + "M"; }
                else if (dlen < 0) {
                    tslen = -dlen; rm += dlen;
                    if (rm < 0) { RDOFF = RDOFF + rm; nc = std::to_string(RDOFF) + "M" + std::to_string(tslen) + "I"; }
                    else { nc += std::to_string(tslen) + "I" + std::to_string(rm) + "M"; }
                } else { nc += std::to_string(dlen) + "D" + std::to_string(rm) + "M"; }
            } else {
                if (dlen == 0) { nc += std::to_string(tslen) + "I" + std::to_string(rm) + "M"; }
                else if (dlen < 0) { rm += dlen; nc += std::to_string(tslen) + "I" + std::to_string(rm) + "M"; }
                else { nc += std::to_string(dlen) + "D" + std::to_string(tslen) + "I" + std::to_string(rm) + "M"; }
            }
            if (mid <= 15) { cs = ov5 + nc + m.suffix().str(); flag = true; }
        }
        // (\d+)D(\d+)M(\d+)([DI])(\d+I)?: combineToCloseToCorrect -- two close deletions (<=15bp gap).
        if (std::regex_search(cs, m, DIG_D_M_DI_I)) {
            int g1 = std::stoi(m[1]), g2 = std::stoi(m[2]), g3 = std::stoi(m[3]); char op = m[4].str()[0];
            if (g2 <= 15) {
                int dlen = g1 + g2, ilen = g2;
                if (op == 'I') ilen += g3;
                else if (op == 'D') { dlen += g3; if (m[5].matched) { std::string istr = m[5].str(); ilen += std::stoi(istr.substr(0, istr.size() - 1)); } }
                cs = m.prefix().str() + std::to_string(dlen) + "D" + std::to_string(ilen) + "I" + m.suffix().str();
                flag = true;
            }
        }
        // (\D)(\d+)I(\d+)M(\d+)([DI])(\d+I)?: combineToCloseToOne -- insertion + short match + indel.
        if (std::regex_search(cs, m, NOTDIG_I_M_DI_I)) {
            std::string lead = m[1].str();
            if (lead != "D" && lead != "H") {
                int g2 = std::stoi(m[2]), g3 = std::stoi(m[3]), g4 = std::stoi(m[4]); char op = m[5].str()[0];
                if (g3 <= 15) {
                    int dlen = g3, ilen = g2 + g3;
                    if (op == 'I') ilen += g4;
                    else if (op == 'D') { dlen += g4; if (m[6].matched) { std::string istr = m[6].str(); ilen += std::stoi(istr.substr(0, istr.size() - 1)); } }
                    std::smatch m2;   // replaceFirst via DIG_I_DIG_M_DIG_DI_DIGI (leading char preserved)
                    if (std::regex_search(cs, m2, DIG_I_M_DI_I))
                        cs = m2.prefix().str() + std::to_string(dlen) + "D" + std::to_string(ilen) + "I" + m2.suffix().str();
                    flag = true;
                }
            }
        }
        // (\d+)D(\d+)D / (\d+)I(\d+)I: merge two adjacent deletions / insertions.
        if (std::regex_search(cs, m, DIG_D_D)) {
            cs = m.prefix().str() + std::to_string(std::stoi(m[1]) + std::stoi(m[2])) + "D" + m.suffix().str();
            flag = true;
        }
        if (std::regex_search(cs, m, DIG_I_I)) {
            cs = m.prefix().str() + std::to_string(std::stoi(m[1]) + std::stoi(m[2])) + "I" + m.suffix().str();
            flag = true;
        }
    }
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
    // complement maps uniquely near the read. CigarModifier.java 83-130 selects the branch by which
    // 2+-digit soft-clip is PRESENT (sc5 = ^(\d\d+)S, sc3 = (\d\d+)S$, i.e. clip length >= 10), with
    // the `>= SEED_2` removal test applied INSIDE the chosen branch. The 5' branch has priority and is
    // an `else if` chain: a 2-digit-but-sub-SEED_2 leading clip (length 10 or 11) therefore claims the
    // 5' branch, does nothing, and suppresses any 3' clip removal. Folding `>= SEED_2` into the branch
    // condition (as before) wrongly let such a read fall through and strip its trailing clip, which in
    // turn hid it from the `^\d\dS.*\d\dS$` chimeric read filter and inflated coverage.
    if (!cfg.chimeric) {
        if (cig.front().second == 'S' && cig.front().first >= 10) {
            int el = cig.front().first;
            if (el >= Reference::SEED_2) {
                int sp = ref.seedUnique(reverseComplement(seq.substr(0, el)).substr(0, Reference::SEED_2));
                if (sp > 0 && std::abs(position - sp) < 2 * maxReadLength) {
                    cig.erase(cig.begin()); seq = seq.substr(el); qual.erase(qual.begin(), qual.begin() + el);
                }
            }
        } else if (cig.back().second == 'S' && cig.back().first >= 10) {
            int el = cig.back().first;
            if (el >= Reference::SEED_2) {
                std::string rc = reverseComplement(seq.substr(seq.size() - el, el));
                int sp = ref.seedUnique(rc.substr(rc.size() - Reference::SEED_2, Reference::SEED_2));
                if (sp > 0 && std::abs(position - sp) < 2 * maxReadLength) {
                    cig.pop_back(); seq = seq.substr(0, seq.size() - el); qual.resize(qual.size() - el);
                }
            }
        }
    }

    // Full port of CigarModifier's `while (flag && indel > 0)` loop: soft-clip normalization of
    // indels wedged next to clips/short matches at read ends, plus collapse of close indel clusters
    // (I+M+D, D+M+D, three indels, adjacent D/D and I/I) into a single D/I complex. Runs on a string
    // CIGAR only when the read carries an indel; the read seq/qual are untouched by the loop.
    {
        int indelLen = 0;
        for (auto& e : cig) if (e.second == 'I' || e.second == 'D') indelLen += e.first;
        if (indelLen > 0) {
            std::string cs = cigToStr(cig);
            indelCollapseLoop(cs, position, seq, ref);
            cig = strToCig(cs);
        }
    }

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
