#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <getopt.h>

#include "config.hpp"
#include "region.hpp"
#include "reference.hpp"
#include "cigar_parser.hpp"
#include "tovars.hpp"
#include "printer.hpp"

using namespace vardict;

static void usage() {
    std::fprintf(stderr,
      "vardictcpp - C++ port of VarDict (simple/pileup counting core)\n"
      "Usage: vardictcpp -G ref.fa -b in.bam -N sample [-R chr:s-e | BED] [options]\n"
      "  -G FILE   indexed reference fasta (required)\n"
      "  -b FILE   indexed BAM (required)\n"
      "  -N STR    sample name\n"
      "  -R STR    region chr:start-end\n"
      "  -c/-S/-E/-g N  BED columns for chr/start/end/gene (1-based; default 1/2/3/4)\n"
      "  -f FLOAT  min allele frequency (default 0.01; 0 keeps all alt positions)\n"
      "  -r INT    min alt reads (default 2)\n"
      "  -q INT    min base quality for hi-qual (default 25)\n"
      "  -O INT    min mapping quality (default 0)\n"
      "  -x INT    region extension bp (default 0)\n"
      "  -p        pileup mode (emit every position)\n"
      "  -t        remove duplicate reads\n"
      "  -h        print header\n"
      "  --chunk INT  split regions longer than INT bp into windows (bounds memory)\n"
      "  (positional: BED file)\n");
}

static std::vector<Region> loadBed(const Config& c) {
    std::vector<Region> regs;
    std::ifstream in(c.bed);
    if (!in) { std::fprintf(stderr, "cannot open BED %s\n", c.bed.c_str()); std::exit(1); }
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::vector<std::string> f;
        std::stringstream ss(line); std::string tok;
        while (std::getline(ss, tok, '\t')) f.push_back(tok);
        int need = std::max({c.colChr, c.colStart, c.colEnd, c.colGene});
        if ((int)f.size() < need) continue;
        Region r;
        r.chr = f[c.colChr - 1];
        int s = std::stoi(f[c.colStart - 1]);
        int e = std::stoi(f[c.colEnd - 1]);
        if (c.zeroBased && s < e) s += 1; // BED half-open -> 1-based
        r.start = s - c.numberNucleotideToExtend;
        r.end   = e + c.numberNucleotideToExtend;
        r.gene = (c.colGene - 1 < (int)f.size()) ? f[c.colGene - 1] : r.chr;
        regs.push_back(r);
    }
    return regs;
}

int main(int argc, char** argv) {
    Config c;
    bool zeroBasedSet = false;

    static struct option longopts[] = {
        {"chunk", required_argument, nullptr, 1000},
        {nullptr, 0, nullptr, 0}
    };
    int ch;
    while ((ch = getopt_long(argc, argv, "G:b:N:R:c:S:E:g:f:r:q:O:x:F:z::pth", longopts, nullptr)) != -1) {
        switch (ch) {
        case 'G': c.ref = optarg; break;
        case 'b': c.bam = optarg; break;
        case 'N': c.sample = optarg; break;
        case 'R': c.region = optarg; break;
        case 'c': c.colChr = std::atoi(optarg); break;
        case 'S': c.colStart = std::atoi(optarg); break;
        case 'E': c.colEnd = std::atoi(optarg); break;
        case 'g': c.colGene = std::atoi(optarg); break;
        case 'f': c.freq = std::atof(optarg); break;
        case 'r': c.minReads = std::atoi(optarg); break;
        case 'q': c.goodq = std::atof(optarg); break;
        case 'O': c.mapqMin = std::atof(optarg); break;
        case 'x': c.numberNucleotideToExtend = std::atoi(optarg); break;
        case 'F': c.samFilterFlag = (int)std::strtol(optarg, nullptr, 0); break;
        case 'z': c.zeroBased = (optarg ? std::atoi(optarg) != 0 : true); zeroBasedSet = true; break;
        case 'p': c.doPileup = true; break;
        case 't': c.removeDuplicates = true; break;
        case 'h': c.printHeader = true; break;
        case 1000: c.chunkSize = std::atoi(optarg); break;
        default: usage(); return 1;
        }
    }
    if (optind < argc) c.bed = argv[optind];
    if (c.ref.empty() || c.bam.empty()) { usage(); return 1; }

    // Build regions.
    std::vector<Region> regions;
    if (!c.region.empty()) {
        Region r;
        auto colon = c.region.find(':');
        r.chr = c.region.substr(0, colon);
        std::string rest = c.region.substr(colon + 1);
        auto dash = rest.find('-');
        r.start = std::stoi(rest.substr(0, dash)) - c.numberNucleotideToExtend;
        r.end   = std::stoi(rest.substr(dash + 1)) + c.numberNucleotideToExtend;
        r.gene = r.chr;
        regions.push_back(r);
    } else if (!c.bed.empty()) {
        // 4-column custom BED defaults to zero-based unless -z given (matches VarDict).
        if (!zeroBasedSet) c.zeroBased = true;
        regions = loadBed(c);
    } else {
        std::fprintf(stderr, "error: need -R or a BED file\n");
        return 1;
    }
    regions = splitLongRegions(regions, c.chunkSize);

    if (c.printHeader) printHeader(stdout);

    Reference ref(c.ref);
    CigarParser parser(c, ref);
    for (const auto& region : regions) {
        // Reference window with a small pad for flank columns.
        ref.load(region.chr, region.start, region.end, 20 + c.numberNucleotideToExtend);
        VariationData vd;
        parser.process(region, vd);
        auto variants = callVariants(c, region, vd, ref);
        for (const auto& v : variants) printVariant(stdout, c, region, v);
        // vd goes out of scope here -> per-region memory is released before the next window.
    }
    return 0;
}
