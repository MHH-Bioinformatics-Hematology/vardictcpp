#pragma once
// Port of modules/CigarModifier.java: reshape a read's CIGAR before counting (mismatch read-ends ->
// soft-clip, chimeric-clip removal). The indel-collapse loop is not ported; reads containing I/D are
// left unchanged. Verified read-by-read against instrumented VarDict (MODCIG diff harness).
#include <vector>
#include <string>
#include <utility>
#include "reference.hpp"
#include "config.hpp"

namespace vardict {

using Cig = std::vector<std::pair<int, char>>; // (length, op)

void modifyCigar(int& position, Cig& cigar, std::string& seq, std::vector<int>& qual,
                 Reference& ref, int maxReadLength, const Config& cfg);

} // namespace vardict
