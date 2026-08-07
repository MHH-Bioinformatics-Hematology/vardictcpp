#pragma once
// Faithful port of com.astrazeneca.vardict.Utils (the Perl-compat string layer) + numeric rounding.
#include <string>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>

namespace vardict {

// substr(string, idx): Perl-style; negative idx counts from the right end.
inline std::string substr(const std::string& s, int idx) {
    if (idx >= 0) return s.substr(std::min((int)s.size(), idx));
    return s.substr(std::max(0, (int)s.size() + idx));
}

// substr(string, begin, len): negative begin counts from right; len<0 means "up to len from the end".
inline std::string substr(const std::string& s, int begin, int len) {
    int n = (int)s.size();
    if (begin < 0) begin = n + begin;
    if (begin < 0) begin = 0;
    if (begin > n) begin = n;
    if (len > 0) {
        int end = std::min(begin + len, n);
        return s.substr(begin, end - begin);
    } else if (len == 0) {
        return "";
    } else {
        int end = n + len;
        if (end < begin) return "";
        return s.substr(begin, end - begin);
    }
}

// charAt with negative-index support; returns (char)-1 (0xFF) when out of range, like Utils.charAt.
inline char charAt(const std::string& s, int index) {
    if (index < 0) {
        int i = (int)s.size() + index;
        if (i < 0) return (char)-1;
        return s[i];
    }
    if (index >= (int)s.size()) return (char)-1;
    return s[index];
}

inline std::string reverseStr(std::string s) { std::reverse(s.begin(), s.end()); return s; }

inline char complementBase(char c) {
    switch (c) {
        case 'A': return 'T'; case 'T': return 'A'; case 'C': return 'G'; case 'G': return 'C';
        case 'a': return 't'; case 't': return 'a'; case 'c': return 'g'; case 'g': return 'c';
        case 'N': return 'N'; default: return c;
    }
}
inline std::string complementStr(std::string s) { for (auto& c : s) c = complementBase(c); return s; }
inline std::string reverseComplement(std::string s) { return complementStr(reverseStr(std::move(s))); }

// roundHalfEven(pattern, value): mimic Java DecimalFormat(pattern).format then parse back to double.
// `pattern` like "0.0000" -> 4 fractional digits, HALF_EVEN (banker's) rounding.
inline double roundHalfEven(const std::string& pattern, double value) {
    int decimals = 0;
    auto dot = pattern.find('.');
    if (dot != std::string::npos) decimals = (int)(pattern.size() - dot - 1);
    // Round the EXACT double value (as Java's DecimalFormat/BigDecimal does), not value*scale.
    // Multiplying by 10^decimals can snap a value that is only epsilon away from a .5 boundary
    // (e.g. 99.0/220.0 == 0.45000000000000001) exactly onto x.5, turning a clean round-up into a
    // false banker's-rounding tie (-> 0.4). glibc printf performs correct round-half-to-even on the
    // exact stored double, matching Java, so format-then-parse reproduces DecimalFormat's result.
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.*f", decimals, value);
    return std::strtod(buf, nullptr);
}

// getRoundedValueToPrint: integer -> "0" pattern; else pattern with trailing zeros stripped.
std::string getRoundedValueToPrint(const std::string& pattern, double value);

// Format a value the way SimpleOutputVariant does: "0" when exactly zero, else fixed precision.
inline std::string fmtOrZero(double v, int decimals) {
    if (v == 0) return "0";
    char buf[64];
    char pat[16];
    std::snprintf(pat, sizeof(pat), "%%.%df", decimals);
    std::snprintf(buf, sizeof(buf), pat, v);
    return buf;
}

} // namespace vardict
