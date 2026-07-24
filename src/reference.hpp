#pragma once
#include <string>
#include <unordered_map>
#include <vector>
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

    // Ensure [start,end] (1-based inclusive) is covered by the loaded window, expanding the loaded
    // seq to the union [min(start,loadedStart), max(end,loadedEnd)] if needed. Mirrors VarDict's
    // ReferenceResource.getReference merging additional windows into referenceSequences (used by the
    // SV path to pull in a far-off deletion breakpoint outside the region window).
    void ensure(int start, int end);

    // 1-based reference base at position p (uppercase); returns 'N' if outside the loaded window.
    char at(int p) const {
        int i = p - loadedStart_;
        if (i < 0 || i >= (int)seq_.size()) return 'N';
        return seq_[i];
    }
    // Whether position p is within the loaded window (mirrors ref.get(p) != null).
    bool has(int p) const { int i = p - loadedStart_; return i >= 0 && i < (int)seq_.size(); }
    int loadedStart() const { return loadedStart_; }
    int loadedEnd() const { return loadedStart_ + (int)seq_.size() - 1; }

    // k-mer -> reference positions (SEED_1=17 and SEED_2=12), built over the loaded window.
    // Mirrors ReferenceResource.addPositionsToSeedSequence; used by findMatch for SV breakpoints.
    const std::vector<int>* seedPositions(const std::string& kmer) const {
        auto it = seed_.find(kmer);
        return it == seed_.end() ? nullptr : &it->second;
    }
    static constexpr int SEED_1 = 17;
    static constexpr int SEED_2 = 12;

private:
    void buildSeed();
    faidx_t* fai_ = nullptr;
    std::string seq_;
    std::string loadedChr_;
    int loadedStart_ = 1;
    std::unordered_map<std::string, std::vector<int>> seed_;
};

} // namespace vardict
