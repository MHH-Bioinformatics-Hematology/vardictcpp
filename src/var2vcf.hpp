#pragma once
#include <string>

namespace vardict {
struct Config;

// Native ports of VarDict's var2vcf_valid.pl (single sample) and var2vcf_paired.pl (tumor|normal),
// converting the fisher-mode variant TSV (as produced by this caller with -F) into VCF. No Perl/R.
std::string var2vcfSingle(const Config& cfg, const std::string& tsv);
std::string var2vcfPaired(const Config& cfg, const std::string& tsv);

} // namespace vardict
