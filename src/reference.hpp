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

    // Ensure [start,end] (1-based inclusive) is covered by the reference, loading a SEPARATE small
    // window [start-pad, end+pad] (matching Java's getReference(modifiedRegion, extension) window,
    // pad = numberNucleotideToExtend + extension) as a DISJOINT segment if [start,end] is not already
    // covered. Mirrors VarDict's ReferenceResource.getReference adding another window to the reference
    // map: Java stores the reference as disjoint [start,end] segments, NOT one contiguous span, so the
    // far-off SV/large-indel breakpoint pulls in only its own small window instead of gap-filling every
    // base between the region and the breakpoint (which could be multiple Mbp -> multi-GB per thread).
    void ensure(int start, int end, int pad);

    // 1-based reference base at position p (uppercase); returns 'N' if outside all loaded windows.
    char at(int p) const {
        int i = p - loadedStart_;
        if (i >= 0 && i < (int)seq_.size()) return seq_[i];   // primary window fast path
        if (extra_.empty()) return 'N';
        return atExtra(p);
    }
    // Whether position p is within any loaded window (mirrors ref.get(p) != null).
    // The trailing SEED_1 bases of a primary window that does NOT reach the contig end are absent
    // from Java's referenceSequences map (ReferenceResource.getReference: siteEnd = exon.length() -
    // SEED_1 unless len == sequenceEnd), so has() returns false there too. Consumers guard at() with
    // has(), so this truncation alone reproduces Java's ref.get(p) == null at the padded window's tail.
    bool has(int p) const {
        int i = p - loadedStart_;
        if (i >= 0 && i < (int)seq_.size() && p <= primaryEffEnd_) return true;   // primary window fast path
        if (extra_.empty()) return false;
        return hasExtra(p);
    }
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
    // A disjoint reference segment loaded by ensure() (1-based, `start` is its first position).
    struct Segment { int start; std::string seq; };
    void buildSeed() const;
    void buildSeedDisjoint() const;   // seed build for the multi-window (post-ensure) case
    void fetchWindow(const std::string& chr, int s, int e);
    // Fetch [s,e] (1-based, s clamped to 1, e clamped to contig by faidx) uppercased; returns the
    // bases and sets `s` to the clamped start. Empty string if outside the contig.
    std::string fetchSeq(const std::string& chr, int& s, int e) const;
    // Slow-path base/coverage lookup in the disjoint extra segments (only reached off the primary
    // window fast path, i.e. positions the SV/large-indel realignment pulled in far from the region).
    char atExtra(int p) const {
        for (const auto& sg : extra_) { int i = p - sg.start; if (i >= 0 && i < (int)sg.seq.size()) return sg.seq[i]; }
        return 'N';
    }
    bool hasExtra(int p) const {
        for (const auto& sg : extra_) { int i = p - sg.start; if (i >= 0 && i < (int)sg.seq.size()) return true; }
        return false;
    }
    // Whether [a,b] lies fully inside one already-loaded extra segment.
    bool inExtra(int a, int b) const {
        for (const auto& sg : extra_) if (a >= sg.start && b <= sg.start + (int)sg.seq.size() - 1) return true;
        return false;
    }
    // Add a freshly-fetched segment [s, s+seq.size()-1], merging it with any overlapping/adjacent
    // existing extra segments so the extra list stays a set of disjoint, position-sorted windows.
    void addSegment(int s, std::string&& seq);
    // Record a genuinely-requested window [a,b] (1-based inclusive); returns true if it added coverage.
    bool addGenuine(int a, int b);
    // Whether the whole span [a,b] lies inside a single genuinely-loaded window.
    bool inGenuine(int a, int b) const {
        for (const auto& iv : genuine_) if (a >= iv.first && b <= iv.second) return true;
        return false;
    }
    faidx_t* fai_ = nullptr;
    std::string seq_;            // primary window loaded by load() (the region + padding)
    std::string loadedChr_;
    int loadedStart_ = 1;        // 1-based first position of seq_ (the primary window)
    int contigLen_ = 0;          // length of loadedChr_ (0 if unknown); used to detect end-of-contig windows
    // Last position of the primary window that Java's referenceSequences map covers. Equal to loadedEnd()
    // when the window reaches the contig end, else loadedEnd() - SEED_1 (Java drops the trailing SEED_1).
    int primaryEffEnd_ = 0;
    // Extra DISJOINT windows pulled in by ensure() for far-off SV/large-indel breakpoints (Java's
    // additional reference-map windows). Kept sorted and non-overlapping. Empty in the common case,
    // so at()/has() take the primary fast path and never scan this list.
    std::vector<Segment> extra_;
    // Genuinely-requested reference windows (1-based inclusive): the primary window plus each ensure()
    // [start,end]. Java's reference map is a set of DISJOINT windows, so the seed index must only cover
    // these genuine windows -- otherwise indexing padding/other segments invents spurious unique k-mers
    // (e.g. a far-off findMatchRev hit that fabricates an <INV>). See ensure()/buildSeed().
    std::vector<std::pair<int,int>> genuine_;
    // Seed index is built lazily on first seedUnique() query and invalidated whenever the loaded window
    // changes. Most regions never hit the SV/large-indel realignment paths that consume it.
    mutable bool seedBuilt_ = false;
    mutable std::unordered_map<uint64_t, int> seed17_, seed12_;  // encoded k-mer -> unique pos (>0) or 0 if repeated
};

} // namespace vardict
