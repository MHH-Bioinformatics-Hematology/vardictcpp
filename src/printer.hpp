#pragma once
#include <cstdio>
#include <string>
#include "config.hpp"
#include "region.hpp"
#include "tovars.hpp"

namespace vardict {

// Emits VarDict simple-mode TSV (mirrors printers/SimpleOutputVariant.java column order).
void printHeader(std::FILE* out);

// Append one variant's TSV line to `out` (used to build a per-region buffer for ordered output).
void appendVariant(std::string& out, const Config& cfg, const Region& region, const Variant& v);

} // namespace vardict
