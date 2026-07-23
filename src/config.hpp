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
    int minReads = 2;         // -r  minimum alt reads to call (conf.minr)
    int minBiasReads = 2;     // -B  minimum reads per strand for the strand-bias flag
    double bias = 0.05;       // conf.bias: per-strand fraction cutoff for the strand-bias flag
    double mapqMin = 0;       // -O  minimum mapping quality (read filter AND isGoodVar threshold)
    double goodq = 22.5;      // -q  base-quality boundary for hi/lo-quality counting + isGoodVar
    int readPosFilter = 5;    // -P  minimum mean read position (isGoodVar)
    double qratio = 1.5;      // -o  minimum hi/lo-quality ratio (isGoodVar)
    double monomerMsiFrequency = 0.25;    // --mfreq
    double nonMonomerMsiFrequency = 0.1;  // --nmfreq
    int vext = 2;             // -X  extension for MNV/complex-variant growth
    int mismatch = 8;         // -m  skip a read whose (NM - indel length) exceeds this (default 8)
    int indelsize = 50;       // -I  indel size / large-indel breakpoint search window
    int minMatch = 0;         // -M
    int SVMINLEN = 1000;      // -L  minimum SV length to spell as <DUP>/<INV>/<DEL>
    static constexpr int EXTENSION = 5000;
    static constexpr int SVMAXLEN = 150000;
    static constexpr int SVFLANK = 50;
    int numberNucleotideToExtend = 0; // -x
    int samFilterFlag = 0x504; // -F
    bool performLocalRealignment = true; // -k (CigarModifier + realignment)
    bool chimeric = false;    // --chimeric
    static constexpr int LOWQUAL = 10;
    bool doPileup = false;    // -p
    bool removeDuplicates = false; // -t
    int threads = 1;          // -th   (region-parallel; core is per-region single-threaded)
    bool printHeader = false; // -h
    bool fisher = false;      // --fisher : add strand-bias Fisher exact p-value + odds-ratio columns
    int chunkSize = 0;        // --chunk : split long regions into windows (bounds memory)
};

} // namespace vardict
