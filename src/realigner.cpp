#include "realigner.hpp"
#include "util.hpp"
#include <algorithm>
#include <vector>
#include <string>

namespace vardict {

// VariationUtils.adjCnt (2-arg form: referenceVar == null). Adds `variant`'s tallies into `varToAdd`.
static void adjCnt(Variation& varToAdd, const Variation& variant) {
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
}

// Insert '&' after the first character (mnt substrings longer than 1 become MNV descriptions).
static std::string decorate(std::string s) {
    if (s.size() > 1) s.insert(1, "&");
    return s;
}

void adjustMNP(VariationData& vd, Reference& /*ref*/, const Config& /*cfg*/, const Region& /*region*/) {
    struct Item { int position; std::string desc; int count; };
    std::vector<Item> tmp;
    for (auto& [pos, m] : vd.mnp)
        for (auto& [desc, cnt] : m)
            tmp.push_back({pos, desc, cnt});
    std::sort(tmp.begin(), tmp.end(), [](const Item& a, const Item& b) {
        if (a.count != b.count) return a.count > b.count;      // count DESC
        if (a.position != b.position) return a.position < b.position; // position ASC
        return a.desc > b.desc;                                 // description DESC
    });

    for (const auto& it : tmp) {
        int position = it.position;
        const std::string& vn = it.desc;
        auto pit = vd.nonInsertionVariants.find(position);
        if (pit == vd.nonInsertionVariants.end()) continue;
        VarMap& varsOnPosition = pit->second;
        auto vrefIt = varsOnPosition.find(vn);
        if (vrefIt == varsOnPosition.end()) continue; // already consumed by indel realignment
        Variation& vref = vrefIt->second;

        std::string mnt = vn;
        auto amp = mnt.find('&');
        if (amp != std::string::npos) mnt.erase(amp, 1); // replaceFirst("&","")

        for (int i = 0; i < (int)mnt.size() - 1; ++i) {
            std::string left = decorate(substr(mnt, 0, i + 1));
            std::string right = decorate(substr(mnt, -((int)mnt.size() - i - 1)));

            // left prefix at the same position
            auto lit = varsOnPosition.find(left);
            if (lit != varsOnPosition.end()) {
                Variation& tref = lit->second;
                if (tref.varsCount > 0 && tref.varsCount < vref.varsCount &&
                    tref.meanPosition / tref.varsCount <= i + 1) {
                    adjCnt(vref, tref);
                    varsOnPosition.erase(left);
                }
            }
            // right suffix at position + i + 1
            auto rposIt = vd.nonInsertionVariants.find(position + i + 1);
            if (rposIt != vd.nonInsertionVariants.end()) {
                auto rit = rposIt->second.find(right);
                if (rit != rposIt->second.end()) {
                    Variation& tref = rit->second;
                    if (tref.varsCount >= 0 && tref.varsCount < vref.varsCount) {
                        adjCnt(vref, tref);
                        vd.refCoverage[position] += tref.varsCount;
                        rposIt->second.erase(right);
                    }
                }
            }
        }
    }
}

} // namespace vardict
