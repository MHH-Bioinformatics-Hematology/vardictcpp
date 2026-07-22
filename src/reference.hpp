#pragma once
#include <string>
#include <htslib/faidx.h>

namespace vardict {

// Indexed FASTA access (mirrors data/ReferenceResource.java, backed by htslib faidx instead of
// htsjdk IndexedFastaSequenceFile). Loads an uppercased window [start,end] (1-based inclusive)
// for a region and answers per-position base queries. Not thread-safe; one per worker thread.
class Reference {
public:
    explicit Reference(const std::string& fastaPath);
    ~Reference();

    // Load reference window for [start,end] (1-based) with `pad` extra bases on each side.
    void load(const std::string& chr, int start, int end, int pad);

    // 1-based reference base at position p (uppercase); returns 'N' if outside the loaded window.
    char at(int p) const {
        int i = p - loadedStart_;
        if (i < 0 || i >= (int)seq_.size()) return 'N';
        return seq_[i];
    }
    int loadedStart() const { return loadedStart_; }
    int loadedEnd() const { return loadedStart_ + (int)seq_.size() - 1; }

private:
    faidx_t* fai_ = nullptr;
    std::string seq_;
    int loadedStart_ = 1;
};

} // namespace vardict
