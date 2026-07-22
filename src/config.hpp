#pragma once
#include <string>

namespace vardict {

// Subset of VarDict's Configuration relevant to the simple/pileup counting core.
// Field names mirror VarDictJava's Configuration.java so the port is auditable.
struct Config {
    std::string bam;          // -b
    std::string ref;          // -G  reference fasta (indexed)
    std::string sample;       // -N
    std::string bed;          // positional BED file (optional)
    std::string region;       // -R  chr:start-end (optional)

    // BED column mapping (1-based, as VarDict), used with a BED file.
    int colChr = 1;           // -c
    int colStart = 2;         // -S
    int colEnd = 3;           // -E
    int colGene = 4;          // -g
    bool zeroBased = false;   // -z (BED default true for 4-col custom; set by caller)

    double freq = 0.01;       // -f  (0 keeps every position with >=1 alt read)
    int minReads = 2;         // -r  minimum alt reads to call
    int minBiasReads = 2;     // -B  minimum reads per strand for the strand-bias flag
    int minMapQ = 0;          // -O  minimum mapping quality
    int minBaseQ = 25;        // -q  phred (VarDict GOODQUALITY=23 default is 25 via -q)
    double qratio = 1.5;      // -o
    int minMatch = 0;         // -M
    int numberNucleotideToExtend = 0; // -x
    int samFilterFlag = 0x504; // -F
    bool doPileup = false;    // -p
    bool removeDuplicates = false; // -t
    int threads = 1;          // -th   (region-parallel; core is per-region single-threaded)
    bool printHeader = false; // -h
    int chunkSize = 0;        // --chunk : split long regions into windows (bounds memory)

    static constexpr int GOODQUALITY = 23; // VarDict Configuration.LOWQUAL boundary is 10; hi-qual base >= this
};

} // namespace vardict
