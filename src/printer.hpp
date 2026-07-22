#pragma once
#include <cstdio>
#include <string>
#include "config.hpp"
#include "region.hpp"
#include "tovars.hpp"

namespace vardict {

// Emits VarDict simple-mode TSV (mirrors printers/SimpleOutputVariant.java column order).
// Streams one line per variant to `out` (no whole-region buffering).
void printHeader(std::FILE* out);
void printVariant(std::FILE* out, const Config& cfg, const Region& region, const Variant& v);

} // namespace vardict
