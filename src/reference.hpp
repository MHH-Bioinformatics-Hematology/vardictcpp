#pragma once
#include <string>
#include <unordered_map>
#include <vector>
#include <cstdint>
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
    // Whether the FASTA index contains a sequence named `name` (used to catch chromosome-naming
    // mismatches such as "chr7" vs "7" before a run silently produces no calls).
    bool hasContig(const std::string& name) const { return fai_ && faidx_has_seq(fai_, name.c_str()); }
    int loadedStart() const { return loadedStart_; }
    int loadedEnd() const { return loadedStart_ + (int)seq_.size() - 1; }

    // Position where `kmer` (length SEED_1=17 or SEED_2=12) occurs UNIQUELY in the loaded window:
    // returns that 1-based position (>0) if the k-mer appears exactly once, 0 if it appears more than
    // once, and -1 if it is absent. Every consumer (findMatch / CigarModifier chimeric-clip) only ever
    // tested `positions.size() == 1` and read positions[0], so this preserves behaviour exactly while
    // storing an encoded-key -> unique-position map instead of a std::string key + std::vector<int> per
    // k-mer occurrence. On a large sparsely-covered -R region (e.g. a whole chromosome) the old index
    // built a heap string for every one of ~150M positions (~20 GB); this is ~4x leaner and byte-identical.
    int seedUnique(const std::string& kmer) const {
        if (!seedBuilt_) buildSeed();
        const auto& m = ((int)kmer.size() == SEED_1) ? seed17_ : seed12_;
        auto it = m.find(encodeKmer(kmer.data(), (int)kmer.size()));
        return it == m.end() ? -1 : it->second;
    }
    // 3 bits/base (A,C,G,T -> 0..3, anything else incl. N -> 4); k <= 21 fits in 63 bits. SEED_1 and
    // SEED_2 are kept in separate maps, so their encodings never collide.
    static uint64_t encodeKmer(const char* s, int k) {
        uint64_t c = 0;
        for (int j = 0; j < k; ++j) {
            char ch = s[j];
            unsigned v = ch=='A'?0u : ch=='C'?1u : ch=='G'?2u : ch=='T'?3u : 4u;
            c = (c << 3) | v;
        }
        return c;
    }
    static constexpr int SEED_1 = 17;
    static constexpr int SEED_2 = 12;

private:
    void buildSeed() const;
    void fetchWindow(const std::string& chr, int s, int e);
    // Record a genuinely-requested window [a,b] (1-based inclusive); returns true if it added coverage.
    bool addGenuine(int a, int b);
    // Whether the whole span [a,b] lies inside a single genuinely-loaded window.
    bool inGenuine(int a, int b) const {
        for (const auto& iv : genuine_) if (a >= iv.first && b <= iv.second) return true;
        return false;
    }
    faidx_t* fai_ = nullptr;
    std::string seq_;
    std::string loadedChr_;
    int loadedStart_ = 1;
    // Genuinely-requested reference windows (1-based inclusive). ensure() gap-fills seq_ contiguously
    // for O(1) base lookups, but Java's reference map is a set of DISJOINT windows, so the seed index
    // must only cover these genuine windows -- otherwise a multi-Mbp gap-fill invents spurious unique
    // k-mers (e.g. a far-off findMatchRev hit that fabricates an <INV>). See ensure()/buildSeed().
    std::vector<std::pair<int,int>> genuine_;
    // Seed index is built lazily on first seedUnique() query and invalidated whenever the loaded window
    // changes. Most regions never hit the SV/large-indel realignment paths that consume it.
    mutable bool seedBuilt_ = false;
    mutable std::unordered_map<uint64_t, int> seed17_, seed12_;  // encoded k-mer -> unique pos (>0) or 0 if repeated
};

} // namespace vardict
