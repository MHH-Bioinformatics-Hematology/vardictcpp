#include "fisher.hpp"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <limits>
#include <algorithm>

namespace vardict {

// ---- Java Utils.roundHalfEven(pattern, value): DecimalFormat(HALF_EVEN).format then parse back.
// snprintf("%.*f") rounds half-to-even under the default FE_TONEAREST mode, matching DecimalFormat.
static double roundHalfEven(int decimals, double value) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.*f", decimals, value);
    return std::atof(buf);
}

// ============================================================================
// commons-math3 SaddlePointExpansion (only the pieces HypergeometricDistribution uses).
// ============================================================================
static const double HALF_LOG_2_PI = 0.5 * std::log(2.0 * M_PI);
static const double TWO_PI = 2.0 * M_PI;

static const double EXACT_STIRLING_ERRORS[] = {
    0.0,                             /* 0.0 */
    0.1534264097200273452913848,     /* 0.5 */
    0.0810614667953272582196702,     /* 1.0 */
    0.0548141210519176538961390,     /* 1.5 */
    0.0413406959554092940938221,     /* 2.0 */
    0.03316287351993628748511048,    /* 2.5 */
    0.02767792568499833914878929,    /* 3.0 */
    0.02374616365629749597132920,    /* 3.5 */
    0.02079067210376509311152277,    /* 4.0 */
    0.01848845053267318523077934,    /* 4.5 */
    0.01664469118982119216319487,    /* 5.0 */
    0.01513497322191737887351255,    /* 5.5 */
    0.01387612882307074799874573,    /* 6.0 */
    0.01281046524292022692424986,    /* 6.5 */
    0.01189670994589177009505572,    /* 7.0 */
    0.01110455975820691732662991,    /* 7.5 */
    0.010411265261972096497478567,   /* 8.0 */
    0.009799416126158803298389475,   /* 8.5 */
    0.009255462182712732917728637,   /* 9.0 */
    0.008768700134139385462952823,   /* 9.5 */
    0.008330563433362871256469318,   /* 10.0 */
    0.007934114564314020547248100,   /* 10.5 */
    0.007573675487951840794972024,   /* 11.0 */
    0.007244554301320383179543912,   /* 11.5 */
    0.006942840107209529865664152,   /* 12.0 */
    0.006665247032707682442354394,   /* 12.5 */
    0.006408994188004207068439631,   /* 13.0 */
    0.006171712263039457647532867,   /* 13.5 */
    0.005951370112758847735624416,   /* 14.0 */
    0.005746216513010115682023589,   /* 14.5 */
    0.005554733551962801371038690    /* 15.0 */
};

static double getStirlingError(double z) {
    double ret;
    if (z < 15.0) {
        double z2 = 2.0 * z;
        if (std::floor(z2) == z2) {
            ret = EXACT_STIRLING_ERRORS[(int) z2];
        } else {
            // z is always integral where called from logBinomialProbability, so this branch
            // (which would need Gamma.logGamma) is never taken here.
            ret = std::lgamma(z + 1.0) - (z + 0.5) * std::log(z) + z - HALF_LOG_2_PI;
        }
    } else {
        double z2 = z * z;
        ret = (0.083333333333333333333 -
                (0.00277777777777777777778 -
                        (0.00079365079365079365079365 -
                                (0.000595238095238095238095238 -
                                        0.0008417508417508417508417508 /
                                        z2) / z2) / z2) / z2) / z;
    }
    return ret;
}

static double getDeviancePart(double x, double mu) {
    double ret;
    if (std::fabs(x - mu) < 0.1 * (x + mu)) {
        double d = x - mu;
        double v = d / (x + mu);
        double s1 = v * d;
        double s = std::numeric_limits<double>::quiet_NaN();
        double ej = 2.0 * x * v;
        v *= v;
        int j = 1;
        while (s1 != s) {
            s = s1;
            ej *= v;
            s1 = s + ej / ((j * 2) + 1);
            ++j;
        }
        ret = s1;
    } else {
        ret = x * std::log(x / mu) + mu - x;
    }
    return ret;
}

static double logBinomialProbability(int x, int n, double p, double q) {
    double ret;
    if (x == 0) {
        if (p < 0.1) {
            ret = -getDeviancePart(n, n * q) - n * p;
        } else {
            ret = n * std::log(q);
        }
    } else if (x == n) {
        if (q < 0.1) {
            ret = -getDeviancePart(n, n * p) - n * q;
        } else {
            ret = n * std::log(p);
        }
    } else {
        ret = getStirlingError(n) - getStirlingError(x) -
              getStirlingError(n - x) - getDeviancePart(x, n * p) -
              getDeviancePart(n - x, n * q);
        double f = (TWO_PI * x * (n - x)) / n;
        ret = -0.5 * std::log(f) + ret;
    }
    return ret;
}

// ============================================================================
// commons-math3 HypergeometricDistribution(populationSize, numberOfSuccesses, sampleSize).
// ============================================================================
struct Hypergeometric {
    int populationSize, numberOfSuccesses, sampleSize;
    int lower, upper;
    Hypergeometric(int pop, int succ, int samp)
        : populationSize(pop), numberOfSuccesses(succ), sampleSize(samp) {
        lower = std::max(0, sampleSize - (populationSize - numberOfSuccesses));
        upper = std::min(numberOfSuccesses, sampleSize);
    }
    double logProbability(int xx) const {
        if (xx < lower || xx > upper) return -std::numeric_limits<double>::infinity();
        double p = (double) sampleSize / (double) populationSize;
        double q = (double) (populationSize - sampleSize) / (double) populationSize;
        double p1 = logBinomialProbability(xx, numberOfSuccesses, p, q);
        double p2 = logBinomialProbability(sampleSize - xx, populationSize - numberOfSuccesses, p, q);
        return p1 + p2;
    }
    double probability(int xx) const {
        double lp = logProbability(xx);
        return lp == -std::numeric_limits<double>::infinity() ? 0.0 : std::exp(lp);
    }
    double innerCumulativeProbability(int x0, int x1, int dx) const {
        double ret = probability(x0);
        while (x0 != x1) {
            x0 += dx;
            ret += probability(x0);
        }
        return ret;
    }
    double cumulativeProbability(int xx) const {
        if (xx < lower) return 0.0;
        if (xx >= upper) return 1.0;
        return innerCumulativeProbability(lower, xx, 1);
    }
    double upperCumulativeProbability(int xx) const {
        if (xx <= lower) return 1.0;
        if (xx > upper) return 0.0;
        return innerCumulativeProbability(upper, xx, -1);
    }
};

// ============================================================================
// UnirootZeroIn.zeroinC (Brent) - direct port.
// ============================================================================
static double zeroinC(double ax, double bx, const std::function<double(double)>& f, double tol) {
    double a, b, c;
    double fa, fb, fc;
    double EPSILON = std::numeric_limits<double>::epsilon(); // Math.ulp(1.0)
    a = ax; b = bx;
    fa = f(a); fb = f(b);
    c = a; fc = fa;
    for (;;) {
        double prev_step = b - a;
        double tol_act;
        double p, q, new_step;
        if (std::fabs(fc) < std::fabs(fb)) {
            a = b; b = c; c = a;
            fa = fb; fb = fc; fc = fa;
        }
        tol_act = 2 * EPSILON * std::fabs(b) + tol / 2.0;
        new_step = (c - b) / 2.0;
        if (std::fabs(new_step) <= tol_act || fb == 0.0) return b;
        if (std::fabs(prev_step) >= tol_act && std::fabs(fa) > std::fabs(fb)) {
            double t1, cb, t2;
            cb = c - b;
            if (a == c) {
                t1 = fb / fa;
                p = cb * t1;
                q = 1.0 - t1;
            } else {
                q = fa / fc;
                t1 = fb / fc;
                t2 = fb / fa;
                p = t2 * (cb * q * (q - t1) - (b - a) * (t1 - 1.0));
                q = (q - 1.0) * (t1 - 1.0) * (t2 - 1.0);
            }
            if (p > 0.0) q = -q;
            else p = -p;
            if (p < (0.75 * cb * q - std::fabs(tol_act * q) / 2.0) && p < std::fabs(prev_step * q / 2.0)) {
                new_step = p / q;
            }
        }
        if (std::fabs(new_step) < tol_act) {
            if (new_step > 0.0) new_step = tol_act;
            else new_step = -tol_act;
        }
        a = b; fa = fb;
        b += new_step;
        fb = f(b);
        if ((fb > 0 && fc > 0) || (fb < 0 && fc < 0)) {
            c = a; fc = fa;
        }
    }
}

// ============================================================================
// FisherExact
// ============================================================================
namespace {
    // support is the contiguous range [lo..hi]; index i corresponds to value lo+i.
    // dnhyper(ncp): the non-central hypergeometric density on the support, normalized.
    std::vector<double> dnhyper(const std::vector<double>& logdc, int lo, int hi, double ncp) {
        int sz = hi - lo + 1;
        std::vector<double> result(sz);
        double logncp = std::log(ncp);
        for (int i = 0; i < sz; ++i) result[i] = logdc[i] + logncp * (lo + i);
        double maxResult = result[0];
        for (int i = 1; i < sz; ++i) if (result[i] > maxResult) maxResult = result[i];
        double sum = 0.0;
        for (int i = 0; i < sz; ++i) { result[i] = std::exp(result[i] - maxResult); sum += result[i]; }
        for (int i = 0; i < sz; ++i) result[i] /= sum;
        return result;
    }
    double mnhyper(const std::vector<double>& logdc, int lo, int hi, double ncp) {
        if (ncp == 0) return (double) lo;
        if (std::isinf(ncp)) return (double) hi;
        std::vector<double> d = dnhyper(logdc, lo, hi, ncp);
        double b = 0.0;
        for (int i = 0; i < (int)d.size(); ++i) b += (double)(lo + i) * d[i];
        return b;
    }
}

FisherExact::FisherExact(int refFwd, int refRev, int altFwd, int altRev) {
    m = refFwd + refRev;
    n = altFwd + altRev;
    k = refFwd + altFwd;
    x = refFwd;
    lo = std::max(0, k - n);
    hi = std::min(k, m);

    // logdcDhyper: density of central hypergeometric on the support, rounded to 7 decimals.
    int sz = hi - lo + 1;
    logdc.resize(sz);
    for (int i = 0; i < sz; ++i) {
        int element = lo + i;
        if (m + n == 0) { logdc[i] = 0.0; continue; }
        Hypergeometric dhyper(m + n, m, k);
        double value = dhyper.logProbability(element);
        if (std::isnan(value)) value = 0.0;
        logdc[i] = roundHalfEven(7, value);
    }

    // calculatePValue: two-sided.
    double relErr = 1 + 1E-7;
    std::vector<double> d = dnhyper(logdc, lo, hi, 1.0);
    double sum = 0.0;
    double pivot = d[x - lo] * relErr;
    for (double el : d) if (el <= pivot) sum += el;
    pvalueTwoSided = sum;
}

// round_as_r: round(value * 1e5) half-even / 1e5, snapping to 0/1.
static double round_as_r(double value) {
    static const double RESULT_ROUND_R = 1E5;
    value = roundHalfEven(0, value * RESULT_ROUND_R);
    value = value / RESULT_ROUND_R;
    value = value == 0.0 ? 0 : (value == 1.0 ? 1 : value);
    return value;
}

double FisherExact::getPValue() const {
    return round_as_r(pvalueTwoSided);
}

// One-sided p-values summed over the NORMALIZED central hypergeometric density (dnhyper), the same
// distribution the two-sided uses. The raw Hypergeometric struct probabilities are unnormalized
// (commons-math binomial-product form), so summing them directly under-counts the tail; normalizing
// first recovers the exact hypergeometric CDF that VarDictJava's pnhyper reports.
double FisherExact::getPValueLess() const {
    if (m + n == 0) return round_as_r(1.0);
    std::vector<double> d = dnhyper(logdc, lo, hi, 1.0);
    double sum = 0.0;
    for (int i = 0; i < (int)d.size(); ++i) if (lo + i <= x) sum += d[i];
    return round_as_r(sum);
}
double FisherExact::getPValueGreater() const {
    if (m + n == 0) return round_as_r(1.0);
    std::vector<double> d = dnhyper(logdc, lo, hi, 1.0);
    double sum = 0.0;
    for (int i = 0; i < (int)d.size(); ++i) if (lo + i >= x) sum += d[i];
    return round_as_r(sum);
}

// mle(x): conditional MLE for the odds ratio.
static double mle(const std::vector<double>& logdc, int lo, int hi, int x) {
    double eps = std::numeric_limits<double>::epsilon(); // Math.ulp(1.0)
    if (x == lo) return 0.0;
    if (x == hi) return std::numeric_limits<double>::infinity();
    double mu = mnhyper(logdc, lo, hi, 1.0);
    double root;
    double tol = std::pow(eps, 0.25);
    if (mu > x) {
        auto f = [&](double t) { return mnhyper(logdc, lo, hi, t) - x; };
        root = zeroinC(0, 1, f, tol);
    } else if (mu < x) {
        auto f = [&](double t) { return mnhyper(logdc, lo, hi, 1.0 / t) - x; };
        root = 1.0 / zeroinC(eps, 1, f, tol);
    } else {
        root = 1.0;
    }
    return root;
}

// String.valueOf(double) for round_as_r outputs (N/1e5 in normal range): shortest decimal,
// which for these values equals the 5-decimal form with trailing zeros stripped.
static std::string formatOddRatioValue(double v) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.5f", v);
    std::string s(buf);
    // strip trailing zeros (keep the decimal digits Java's shortest repr would keep)
    size_t last = s.find_last_not_of('0');
    if (s[last] == '.') last++; // keep one digit after dot -> "x.0" (Java never emits a bare "x.")
    s.erase(last + 1);
    return s;
}

std::string FisherExact::getOddRatio() const {
    double oddRatio = mle(logdc, lo, hi, x);
    if (std::isinf(oddRatio)) return "Inf";
    // value == Math.round(value): integral -> DecimalFormat("0")
    if (oddRatio == std::floor(oddRatio + 0.5)) {
        char buf[64];
        std::snprintf(buf, sizeof(buf), "%.0f", oddRatio);
        return buf;
    }
    return formatOddRatioValue(round_as_r(oddRatio));
}

} // namespace vardict
