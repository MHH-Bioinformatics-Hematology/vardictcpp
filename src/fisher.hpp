#pragma once
#include <string>
#include <vector>

namespace vardict {

// Port of data/fishertest/FisherExact.java (+ the commons-math3 HypergeometricDistribution and
// SaddlePointExpansion it depends on, and UnirootZeroIn). Two-sided Fisher exact test on the
// 2x2 strand table (refFwd,refRev,altFwd,altRev), matching R's fisher.test as VarDict reimplements it.
class FisherExact {
public:
    FisherExact(int refFwd, int refRev, int altFwd, int altRev);
    double getPValue() const;      // two-sided p-value, rounded as R (5 decimals)
    std::string getOddRatio() const; // conditional MLE odds ratio, formatted as VarDict prints it

private:
    int m, n, k, x;
    int lo, hi;
    double pvalueTwoSided = 0;
    std::vector<double> logdc; // stored on the support [lo..hi]
};

} // namespace vardict
