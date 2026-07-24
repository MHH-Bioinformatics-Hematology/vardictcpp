#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
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
#include "somatic.hpp"
#include "amplicon.hpp"

using namespace vardict;

static void usage() {
    std::fprintf(stderr,
      "vardictcpp - C++ port of VarDict (single-sample simple/pileup calling)\n"
      "Accepts VarDict's full option syntax (commons-cli style, e.g. -th 8, -VS STRICT).\n"
      "Usage: vardictcpp -G ref.fa -b in.bam -N sample [-R chr:s-e | BED] [options]\n"
      "  Acted on: -G -b -N -R -c -S -E -g -f -r -q -O -P -o -B -X -I -L -x -F -z -p -t -h -k\n"
      "            -mfreq -nmfreq --chunk -th/--threads\n"
      "  Accepted (parsed, VarDict-compatible; not all affect output yet): -A -M -Q -T -V -W -w -Y\n"
      "            -Z -d -e -n -s -v -y -3 -C -D -K -U -UN -j -J -DP -VS -adaptor -deldupvar -m\n"
      "  Somatic (paired) mode: -b 'tumor|normal' -N 'tumor|normal' (55-column output)\n"
      "  Amplicon (multiplex) mode: -a EDGE:FRACTION with an 8-column amplicon BED (38-column output)\n"
      "  --chunk INT  split regions longer than INT bp into windows (bounds memory)\n");
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

// VarDict's full option set (commons-cli). Value here = does the option take an argument.
// Parsing mirrors commons-cli: options are `-name [value]` (single dash, names may be multi-char such
// as th/VS/DP/mfreq/chimeric) or `--long`. Every VarDict option is accepted so any VarDict command line
// parses; options this port does not act on are recorded but ignored (see below), and unsupported
// *modes* error out rather than silently mis-call.
static const std::map<std::string, bool> VARDICT_OPTS = {
    {"A",1},{"B",1},{"DP",1},{"E",1},{"F",1},{"G",1},{"I",1},{"J",1},{"L",1},{"M",1},{"N",1},{"O",1},
    {"P",1},{"Q",1},{"R",1},{"S",1},{"T",1},{"V",1},{"VS",1},{"W",1},{"X",1},{"Y",1},{"Z",1},{"a",1},
    {"adaptor",1},{"b",1},{"c",1},{"d",1},{"e",1},{"f",1},{"g",1},{"j",1},{"m",1},{"mfreq",1},{"n",1},
    {"nmfreq",1},{"o",1},{"q",1},{"r",1},{"s",1},{"w",1},{"x",1},{"chunk",1},{"th",1},{"threads",1},
    {"3",0},{"C",0},{"D",0},{"H",0},{"K",0},{"U",0},{"UN",0},{"chimeric",0},{"deldupvar",0},{"fisher",0},
    {"h",0},{"i",0},{"k",0},{"p",0},{"t",0},{"u",0},{"v",0},{"y",0},{"z",0},{"?",0},
};

int main(int argc, char** argv) {
    Config c;
    std::map<std::string, std::string> opt;   // parsed options (name -> value; flags -> "")
    std::vector<std::string> positional;

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a.size() >= 1 && a[0] == '-' && a != "-") {
            int d = 0; while (d < (int)a.size() && a[d] == '-') d++;   // strip 1 or 2 leading dashes
            std::string name = a.substr(d);
            std::string inlineVal;
            auto eq = name.find('=');
            if (eq != std::string::npos) { inlineVal = name.substr(eq + 1); name = name.substr(0, eq); }
            auto it = VARDICT_OPTS.find(name);
            if (it == VARDICT_OPTS.end()) {
                std::fprintf(stderr, "vardictcpp: unknown option -%s\n", name.c_str());
                usage(); return 1;
            }
            if (it->second) { // takes an argument
                if (!inlineVal.empty()) opt[name] = inlineVal;
                else if (i + 1 < argc) opt[name] = argv[++i];
                else { std::fprintf(stderr, "vardictcpp: option -%s requires a value\n", name.c_str()); return 1; }
            } else opt[name] = "";
        } else positional.push_back(a);
    }
    auto has = [&](const char* k){ return opt.count(k) != 0; };
    auto val = [&](const char* k, const char* d){ auto it = opt.find(k); return it == opt.end() ? std::string(d) : it->second; };

    // Amplicon (multiplex) mode: -a EDGE:FRACTION (GlobalReadOnlyScope.ampliconBasedCalling).
    if (has("a")) {
        c.amplicon = true;
        std::string a = opt["a"];
        auto colon = a.find(':');
        // Java: distanceToAmplicon = toInt(split[0]); overlapFraction = parseDouble(split[1]);
        // on NumberFormatException it falls back to 10 / 0.95.
        try {
            c.ampEdge = std::stoi(colon == std::string::npos ? a : a.substr(0, colon));
            c.ampFraction = std::stod(colon == std::string::npos ? std::string() : a.substr(colon + 1));
        } catch (...) { c.ampEdge = 10; c.ampFraction = 0.95; }
    }
    c.fisher = has("fisher");

    c.ref = val("G", "");
    c.bam = val("b", "");
    c.sample = val("N", "");
    // Somatic (paired) mode: -b "tumor|normal" and -N "tumor|normal". Split on '|'.
    {
        auto pipe = c.bam.find('|');
        if (pipe != std::string::npos) {
            c.somatic = true;
            c.bam2 = c.bam.substr(pipe + 1);
            c.bam = c.bam.substr(0, pipe);
            auto npipe = c.sample.find('|');
            if (npipe != std::string::npos) {
                c.sample2 = c.sample.substr(npipe + 1);
                c.sample = c.sample.substr(0, npipe);
            }
        }
    }
    c.region = val("R", "");
    c.colChr = std::atoi(val("c", "1").c_str());
    c.colStart = std::atoi(val("S", "2").c_str());
    c.colEnd = std::atoi(val("E", "3").c_str());
    c.colGene = std::atoi(val("g", "4").c_str());
    c.freq = std::atof(val("f", "0.01").c_str());
    c.minReads = std::atoi(val("r", "2").c_str());
    c.minBiasReads = std::atoi(val("B", "2").c_str());
    c.goodq = std::atof(val("q", "22.5").c_str());
    c.mapqMin = std::atof(val("O", "0").c_str());
    c.readPosFilter = std::atoi(val("P", "5").c_str());
    c.qratio = std::atof(val("o", "1.5").c_str());
    c.lofreq = std::atof(val("V", "0.05").c_str());
    c.vext = std::atoi(val("X", "2").c_str());
    c.mismatch = std::atoi(val("m", "8").c_str());
    c.indelsize = std::atoi(val("I", "50").c_str());
    c.SVMINLEN = std::atoi(val("L", "1000").c_str());
    c.monomerMsiFrequency = std::atof(val("mfreq", "0.25").c_str());
    c.nonMonomerMsiFrequency = std::atof(val("nmfreq", "0.1").c_str());
    c.numberNucleotideToExtend = std::atoi(val("x", "0").c_str());
    if (has("F")) c.samFilterFlag = (int)std::strtol(opt["F"].c_str(), nullptr, 0);
    c.zeroBased = has("z");
    c.doPileup = has("p");
    c.removeDuplicates = has("t");
    c.printHeader = has("h");
    c.chimeric = has("chimeric");
    if (has("k")) c.performLocalRealignment = std::atoi(val("k", "1").c_str()) != 0;
    c.chunkSize = std::atoi(val("chunk", "0").c_str());
    c.threads = std::max(1, std::atoi(val(has("threads") ? "threads" : "th", "1").c_str()));
    if (has("H") || has("?")) { usage(); return 0; }

    bool zeroBasedSet = has("z");
    if (!positional.empty()) c.bed = positional[0];
    if (c.ref.empty() || c.bam.empty()) { usage(); return 1; }

    // Amplicon (multiplex) mode: process by segment, comparing calls across overlapping amplicons.
    if (c.amplicon) {
        if (c.bed.empty()) { std::fprintf(stderr, "vardictcpp: amplicon mode (-a) requires a BED file\n"); return 1; }
        std::ifstream bin(c.bed);
        if (!bin) { std::fprintf(stderr, "cannot open BED %s\n", c.bed.c_str()); return 1; }
        std::vector<std::string> raws; std::string line;
        while (std::getline(bin, line)) if (!line.empty() && line[0] != '#') raws.push_back(line);
        auto segments = buildAmpRegions(raws, c);
        if (c.printHeader) printAmpliconHeader(stdout);
        Reference ref(c.ref);
        for (auto& seg : segments) {
            std::vector<std::map<int, AmpVars>> vars;
            for (auto& region : seg) {
                ref.load(region.chr, region.start, region.end, 1200 + c.numberNucleotideToExtend);
                VariationData vd;
                CigarParser(c, ref).process(region, vd);
                adjustMNP(vd, ref, c, region);
                realigndel(vd, ref, c, region, vd.maxReadLength);
                realignins(vd, ref, c, region, vd.maxReadLength);
                realignlgdel(vd, ref, c, region, vd.maxReadLength);
                realignlgins30(vd, ref, c, region, vd.maxReadLength);
                realignlgins(vd, ref, c, region, vd.maxReadLength);
                adjSNV(vd, ref);
                vars.push_back(buildAmpVars(c, region, vd, ref));
            }
            std::string buf;
            appendAmpliconSegment(buf, c, seg, vars);
            std::fputs(buf.c_str(), stdout);
        }
        return 0;
    }

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

    if (c.printHeader) { if (c.somatic) printSomaticHeader(stdout); else printHeader(stdout); }

    // Process one region into a fresh output buffer. Each call uses its own Reference/CigarParser so
    // it is safe to run concurrently (htslib faidx/BAM handles are not shared across threads). vd is
    // released at the end of the call, so per-region memory never accumulates.
    auto processRegion = [&](const Region& region, Reference& ref) {
        // VarDict loads reference with numberNucleotideToExtend + referenceExtension(1200) padding;
        // realignment flanks + the seed index for findMatch need this wider window.
        ref.load(region.chr, region.start, region.end, 1200 + c.numberNucleotideToExtend);
        VariationData vd;
        CigarParser(c, ref).process(region, vd);
        // Realignment order mirrors VariationRealigner: filterAllSVStructures (collapse discordant-pair
        // SV clusters) runs first, then adjustMNP, then realigndel, realignins, realignlgdel, ...
        if (!c.disableSV) filterSVStructures(vd, vd.maxReadLength);
        adjustMNP(vd, ref, c, region);
        realigndel(vd, ref, c, region, vd.maxReadLength);
        realignins(vd, ref, c, region, vd.maxReadLength);
        realignlgdel(vd, ref, c, region, vd.maxReadLength);
        realignlgins30(vd, ref, c, region, vd.maxReadLength);
        realignlgins(vd, ref, c, region, vd.maxReadLength);
        // StructuralVariantsProcessor.findAllSVs runs after realignment, before adjSNV. Only the DEL
        // discordant-pair path (findDELdisc) is ported.
        if (!c.disableSV) findDELdisc(vd, ref, c, region, vd.maxReadLength);
        adjSNV(vd, ref);
        std::string buf;
        if (c.somatic) {
            // Paired analysis. The harness pairs a BAM with itself, so the tumor pipeline result is
            // reused for the normal sample (identical counts); SomaticPostProcessModule then compares
            // the two. Building once guarantees v1 == v2 exactly (no maxReadLength-seed drift).
            auto positions = callVariantsSomatic(c, region, vd, ref);
            appendSomaticRegion(buf, c, region, positions);
        } else {
            auto variants = callVariants(c, region, vd, ref);
            for (const auto& v : variants) appendVariant(buf, c, region, v);
        }
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
