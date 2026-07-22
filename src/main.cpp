#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <atomic>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <getopt.h>

#include "config.hpp"
#include "region.hpp"
#include "reference.hpp"
#include "cigar_parser.hpp"
#include "realigner.hpp"
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
        {"chunk",   required_argument, nullptr, 1000},
        {"threads", required_argument, nullptr, 1001},
        {"th",      required_argument, nullptr, 1001},
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
        case 1001: c.threads = std::max(1, std::atoi(optarg)); break;
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
        // VarDict only auto-zero-bases the 4-column *custom* format, which applies when -c is NOT set.
        // When -c/-S/-E/-g are given the BED is treated 1-based unless -z is passed explicitly.
        (void)zeroBasedSet;
        regions = loadBed(c);
    } else {
        std::fprintf(stderr, "error: need -R or a BED file\n");
        return 1;
    }
    regions = splitLongRegions(regions, c.chunkSize);

    if (c.printHeader) printHeader(stdout);

    // Process one region into a fresh output buffer. Each call uses its own Reference/CigarParser so
    // it is safe to run concurrently (htslib faidx/BAM handles are not shared across threads). vd is
    // released at the end of the call, so per-region memory never accumulates.
    auto processRegion = [&](const Region& region, Reference& ref) {
        // VarDict loads reference with numberNucleotideToExtend + referenceExtension(1200) padding;
        // realignment flanks + the seed index for findMatch need this wider window.
        ref.load(region.chr, region.start, region.end, 1200 + c.numberNucleotideToExtend);
        VariationData vd;
        CigarParser(c, ref).process(region, vd);
        // Realignment order mirrors VariationRealigner: adjustMNP, then realigndel, realignins,
        // realignlgdel (large-deletion soft-clip breakpoints).
        adjustMNP(vd, ref, c, region);
        realigndel(vd, ref, c, region, vd.maxReadLength);
        realignins(vd, ref, c, region, vd.maxReadLength);
        realignlgdel(vd, ref, c, region, vd.maxReadLength);
        realignlgins30(vd, ref, c, region, vd.maxReadLength);
        auto variants = callVariants(c, region, vd, ref);
        std::string buf;
        for (const auto& v : variants) appendVariant(buf, c, region, v);
        return buf;
    };

    const int nreg = (int)regions.size();
    int nthreads = std::max(1, std::min(c.threads, nreg));
    if (nthreads == 1) {
        Reference ref(c.ref);
        for (const auto& region : regions) std::fputs(processRegion(region, ref).c_str(), stdout);
        return 0;
    }

    // Region-parallel with ordered streaming output (mirrors VarDictJava's AbstractParallelMode:
    // workers steal regions; a single consumer flushes buffers in region order so output is
    // deterministic and memory stays bounded to a small window of completed-but-unprinted regions).
    std::vector<std::string> results(nreg);
    std::vector<char> done(nreg, 0);
    std::atomic<int> nextWork{0};
    std::mutex m;
    std::condition_variable cv;
    int printed = 0;

    auto worker = [&]() {
        Reference ref(c.ref); // per-thread faidx handle
        int i;
        while ((i = nextWork.fetch_add(1)) < nreg) {
            std::string buf = processRegion(regions[i], ref);
            std::lock_guard<std::mutex> lk(m);
            results[i] = std::move(buf);
            done[i] = 1;
            cv.notify_all();
        }
    };

    std::vector<std::thread> pool;
    for (int t = 0; t < nthreads; ++t) pool.emplace_back(worker);

    // Consumer: flush region buffers strictly in order as they complete.
    {
        std::unique_lock<std::mutex> lk(m);
        while (printed < nreg) {
            cv.wait(lk, [&]{ return done[printed]; });
            std::string out = std::move(results[printed]);
            results[printed].clear(); results[printed].shrink_to_fit();
            done[printed] = 0;
            ++printed;
            lk.unlock();
            std::fputs(out.c_str(), stdout);
            lk.lock();
        }
    }
    for (auto& t : pool) t.join();
    return 0;
}
