#include "reference.hpp"
#include <stdexcept>
#include <cctype>
#include <cstdlib>
#include <cstring>

namespace vardict {

Reference::Reference(const std::string& fastaPath) {
    fai_ = fai_load(fastaPath.c_str());
    if (!fai_) throw std::runtime_error("cannot load faidx for " + fastaPath +
                                        " (is it indexed with `samtools faidx`?)");
}

Reference::~Reference() {
    if (fai_) fai_destroy(fai_);
}

void Reference::load(const std::string& chr, int start, int end, int pad) {
    int s = start - pad; if (s < 1) s = 1;
    int e = end + pad;
    int len = 0;
    // faidx_fetch_seq is 0-based, end-inclusive.
    char* raw = faidx_fetch_seq(fai_, chr.c_str(), s - 1, e - 1, &len);
    if (!raw || len < 0) {
        // Region outside contig / not found: leave empty so at() returns 'N'.
        seq_.clear();
        loadedStart_ = s;
        if (raw) free(raw);
        return;
    }
    seq_.assign(raw, len);
    for (auto& c : seq_) c = (char)std::toupper((unsigned char)c);
    loadedStart_ = s;
    free(raw);
}

} // namespace vardict
