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
#include <stdexcept>
#include <set>
#include <memory>

#include "config.hpp"
#include "region.hpp"
#include "reference.hpp"
#include "simd.hpp"

#ifndef VARDICTCPP_VERSION
#define VARDICTCPP_VERSION "1"
#endif
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

// Parse an integer option/field value, rejecting anything that is not a whole number so a typo turns
// into a clear message instead of a silent 0 (std::atoi) or an uncaught std::invalid_argument crash.
// `what` is woven into the error, e.g. "option -f" or "BED line 12 start".
static long parseIntOr(const std::string& v, const std::string& what) {
    try {
        size_t pos = 0;
        long r = std::stol(v, &pos);
        if (pos != v.size()) throw std::invalid_argument(v);
        return r;
    } catch (const std::exception&) {
        throw std::runtime_error(what + " expects an integer, got '" + v + "'");
    }
}

static double parseDoubleOr(const std::string& v, const std::string& what) {
    try {
        size_t pos = 0;
        double r = std::stod(v, &pos);
        if (pos != v.size()) throw std::invalid_argument(v);
        return r;
    } catch (const std::exception&) {
        throw std::runtime_error(what + " expects a number, got '" + v + "'");
    }
}

static std::vector<Region> loadBed(const Config& c) {
    std::vector<Region> regs;
    std::ifstream in(c.bed);
    if (!in) throw std::runtime_error("cannot open BED file '" + c.bed + "'");
    const int need = std::max({c.colChr, c.colStart, c.colEnd, c.colGene});
    std::string line;
    long lineno = 0;
    while (std::getline(in, line)) {
        ++lineno;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line[0] == '#') continue;
        if (line.rfind("track", 0) == 0 || line.rfind("browser", 0) == 0) continue;
        std::vector<std::string> f;
        std::stringstream ss(line); std::string tok;
        while (std::getline(ss, tok, '\t')) f.push_back(tok);
        const std::string where = "BED '" + c.bed + "' line " + std::to_string(lineno);
        if ((int)f.size() < need) {
            throw std::runtime_error(where + ": found " + std::to_string(f.size()) +
                " column(s) but -c/-S/-E/-g need at least " + std::to_string(need) +
                " (the file must be tab-separated; supply a gene/name column or adjust -g)");
        }
        Region r;
        r.chr = f[c.colChr - 1];
        if (r.chr.empty()) throw std::runtime_error(where + ": empty chromosome name in column " + std::to_string(c.colChr));
        int s = (int)parseIntOr(f[c.colStart - 1], where + " start (column " + std::to_string(c.colStart) + ")");
        int e = (int)parseIntOr(f[c.colEnd - 1], where + " end (column " + std::to_string(c.colEnd) + ")");
        if (s > e) throw std::runtime_error(where + ": start " + std::to_string(s) + " is greater than end " + std::to_string(e));
        if (c.zeroBased && s < e) s += 1; // BED half-open -> 1-based
        r.start = s - c.numberNucleotideToExtend;
        r.end   = e + c.numberNucleotideToExtend;
        r.gene = (c.colGene - 1 < (int)f.size()) ? f[c.colGene - 1] : r.chr;
        regs.push_back(r);
    }
    if (regs.empty()) throw std::runtime_error("no usable regions found in BED file '" + c.bed + "'");
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

static int run(int argc, char** argv) {
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--version" || a == "-version") {
            std::printf("vardictcpp %s (SIMD backend: %s)\n", VARDICTCPP_VERSION, simd::backend());
            return 0;
        }
    }
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
    // Validated numeric accessors: a bad user value (e.g. "-f abc") becomes a clear message naming the
    // option instead of silently parsing to 0. Defaults are always well-formed, so they never throw.
    auto ival = [&](const char* k, const char* d){ return (int)parseIntOr(val(k, d), std::string("option -") + k); };
    auto dval = [&](const char* k, const char* d){ return parseDoubleOr(val(k, d), std::string("option -") + k); };

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
    c.colChr = ival("c", "1");
    c.colStart = ival("S", "2");
    c.colEnd = ival("E", "3");
    c.colGene = ival("g", "4");
    c.freq = dval("f", "0.01");
    c.minReads = ival("r", "2");
    c.minBiasReads = ival("B", "2");
    c.goodq = dval("q", "22.5");
    c.mapqMin = dval("O", "0");
    c.readPosFilter = ival("P", "5");
    c.qratio = dval("o", "1.5");
    c.lofreq = dval("V", "0.05");
    c.vext = ival("X", "2");
    c.mismatch = ival("m", "8");
    c.indelsize = ival("I", "50");
    c.SVMINLEN = ival("L", "1000");
    c.monomerMsiFrequency = dval("mfreq", "0.25");
    c.nonMonomerMsiFrequency = dval("nmfreq", "0.1");
    c.numberNucleotideToExtend = ival("x", "0");
    if (has("F")) c.samFilterFlag = (int)std::strtol(opt["F"].c_str(), nullptr, 0);
    c.zeroBased = has("z");
    c.doPileup = has("p");
    c.removeDuplicates = has("t");
    c.printHeader = has("h");
    c.chimeric = has("chimeric");
    if (has("k")) c.performLocalRealignment = std::atoi(val("k", "1").c_str()) != 0;
    c.chunkSize = ival("chunk", "0");
    c.threads = std::max(1, ival(has("threads") ? "threads" : "th", "1"));
    if (has("H") || has("?")) { usage(); return 0; }

    bool zeroBasedSet = has("z");
    if (!positional.empty()) c.bed = positional[0];
    if (c.ref.empty()) { std::fprintf(stderr, "vardictcpp: a reference FASTA is required (-G ref.fa)\n\n"); usage(); return 1; }
    if (c.bam.empty()) { std::fprintf(stderr, "vardictcpp: an input BAM is required (-b in.bam)\n\n"); usage(); return 1; }
    if (c.sample.empty()) { std::fprintf(stderr, "vardictcpp: a sample name is required (-N sample)\n\n"); usage(); return 1; }
    if (!c.region.empty() && !c.bed.empty())
        std::fprintf(stderr, "vardictcpp: warning: both a region (-R) and a BED were given; using -R and ignoring the BED\n");

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
        BamReader bam(c.bam);
        for (auto& seg : segments) {
            std::vector<std::map<int, AmpVars>> vars;
            for (auto& region : seg) {
                ref.load(region.chr, region.start, region.end, 1200 + c.numberNucleotideToExtend);
                VariationData vd;
                CigarParser(c, ref, bam).process(region, vd);
                adjustMNP(vd, ref, c, region);
                auto reload = [&](int ms, int me) {
                    Region rr; rr.chr = region.chr; rr.start = ms - 200; rr.end = me + 200; rr.gene = region.gene;
                    CigarParser(c, ref, bam).process(rr, vd, /*reloadMode=*/true);
                };
                realigndel(vd, ref, c, region, vd.maxReadLength, {&bam});
                realignins(vd, ref, c, region, vd.maxReadLength);
                realignlgdel(vd, ref, c, region, vd.maxReadLength, reload, {&bam});
                realignlgins30(vd, ref, c, region, vd.maxReadLength, {&bam});
                realignlgins(vd, ref, c, region, vd.maxReadLength, reload, {&bam});
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
        // Expect chr:start-end (commas allowed, e.g. chr7:55,019,017-55,211,628). A bare "chr7" or any
        // other shape is a common mistake, so report it clearly instead of crashing in the parser.
        auto colon = c.region.find(':');
        auto dash = (colon == std::string::npos) ? std::string::npos : c.region.find('-', colon + 1);
        if (colon == std::string::npos || dash == std::string::npos || colon == 0)
            throw std::runtime_error("invalid region '" + c.region +
                "' for -R; expected chr:start-end, e.g. chr7:1-159138663");
        Region r;
        r.chr = c.region.substr(0, colon);
        std::string ss = c.region.substr(colon + 1, dash - colon - 1);
        std::string es = c.region.substr(dash + 1);
        ss.erase(std::remove(ss.begin(), ss.end(), ','), ss.end());
        es.erase(std::remove(es.begin(), es.end(), ','), es.end());
        int s = (int)parseIntOr(ss, "region start in -R '" + c.region + "'");
        int e = (int)parseIntOr(es, "region end in -R '" + c.region + "'");
        if (s > e) throw std::runtime_error("region start " + std::to_string(s) +
            " is greater than end " + std::to_string(e) + " in -R '" + c.region + "'");
        r.start = s - c.numberNucleotideToExtend;
        r.end   = e + c.numberNucleotideToExtend;
        r.gene = r.chr;
        regions.push_back(r);
    } else if (!c.bed.empty()) {
        // VarDict only auto-zero-bases the 4-column *custom* format, which applies when -c is NOT set.
        // When -c/-S/-E/-g are given the BED is treated 1-based unless -z is passed explicitly.
        (void)zeroBasedSet;
        regions = loadBed(c);
    } else {
        std::fprintf(stderr, "vardictcpp: no target given; pass a region (-R chr:start-end) or a BED file\n\n");
        usage();
        return 1;
    }
    regions = splitLongRegions(regions, c.chunkSize);

    // Preflight: open the reference and BAM once (surfacing missing-file / missing-index errors up
    // front with a clear message) and check that the requested chromosomes actually exist in both.
    // Mismatched naming (e.g. "chr7" vs "7") is the usual reason a run silently produces no calls.
    {
        Reference refChk(c.ref);
        BamReader bamChk(c.bam);
        std::set<std::string> chrs;
        for (const auto& r : regions) chrs.insert(r.chr);
        size_t missingRef = 0;
        for (const auto& chr : chrs) {
            bool inRef = refChk.hasContig(chr);
            bool inBam = bamChk.hasContig(chr);
            if (!inRef) {
                ++missingRef;
                std::fprintf(stderr, "vardictcpp: warning: chromosome '%s' is not in the reference '%s'; it will be skipped\n",
                             chr.c_str(), c.ref.c_str());
            } else if (!inBam) {
                std::fprintf(stderr, "vardictcpp: warning: chromosome '%s' is not in the BAM '%s'; no reads there\n",
                             chr.c_str(), c.bam.c_str());
            }
        }
        if (missingRef == chrs.size())
            throw std::runtime_error("none of the requested chromosomes were found in the reference '" + c.ref +
                "'; check that the chromosome naming matches (e.g. 'chr7' vs '7')");
    }

    if (c.printHeader) { if (c.somatic) printSomaticHeader(stdout); else printHeader(stdout); }

    // Process one region into a fresh output buffer. Each call uses its own Reference/CigarParser so
    // it is safe to run concurrently (htslib faidx/BAM handles are not shared across threads). vd is
    // released at the end of the call, so per-region memory never accumulates.
    // Run the full counting + realignment pipeline for one BAM over a region and return its per-position
    // variation data. VarDict loads reference with numberNucleotideToExtend + referenceExtension(1200)
    // padding; realignment flanks + the seed index for findMatch need this wider window.
    auto runPipeline = [&](Reference& ref, BamReader& b, const Region& region) {
        ref.load(region.chr, region.start, region.end, 1200 + c.numberNucleotideToExtend);
        VariationData vd;
        CigarParser(c, ref, b).process(region, vd);
        // Realignment order mirrors VariationRealigner: filterAllSVStructures (collapse discordant-pair
        // SV clusters) runs first, then adjustMNP, then realigndel, realignins, realignlgdel, ...
        if (!c.disableSV) filterSVStructures(vd, vd.maxReadLength);
        adjustMNP(vd, ref, c, region);
        // reload(ms,me) re-reads coverage over [ms-200, me+200] into vd (reloadMode: no SV clusters).
        auto reload = [&](int ms, int me) {
            Region rr; rr.chr = region.chr; rr.start = ms - 200; rr.end = me + 200; rr.gene = region.gene;
            CigarParser(c, ref, b).process(rr, vd, /*reloadMode=*/true);
        };
        realigndel(vd, ref, c, region, vd.maxReadLength, {&b});
        realignins(vd, ref, c, region, vd.maxReadLength);
        realignlgdel(vd, ref, c, region, vd.maxReadLength, reload, {&b});
        realignlgins30(vd, ref, c, region, vd.maxReadLength, {&b});
        realignlgins(vd, ref, c, region, vd.maxReadLength, reload, {&b});
        // StructuralVariantsProcessor.findAllSVs runs after realignment, before adjSNV. Ported paths,
        // in Java order: findINV (pair-assisted <INV>), findsv (split-read <INV>), findDELdisc (<DEL>).
        if (!c.disableSV) findDEL(vd, ref, c, region, vd.maxReadLength, reload);
        if (!c.disableSV) findINV(vd, ref, c, region, vd.maxReadLength, reload, {&b});
        if (!c.disableSV) findsv(vd, ref, c, region, vd.maxReadLength);
        if (!c.disableSV) findDELdisc(vd, ref, c, region, vd.maxReadLength);
        adjSNV(vd, ref);
        return vd;
    };

    // combineAnalysis merged pipeline: CigarParser runs over BOTH BAMs into ONE VariationData (they
    // accumulate), then the same realign/SV/adjSNV sequence as runPipeline. Mirrors SomaticMode.pipeline
    // invoked with bam = bam1:bam2 inside SomaticPostProcessModule.combineAnalysis.
    auto runMergedPipeline = [&](Reference& ref, BamReader& b1, BamReader& b2, const Region& region) {
        ref.load(region.chr, region.start, region.end, 1200 + c.numberNucleotideToExtend);
        VariationData vd;
        CigarParser(c, ref, b1).process(region, vd);
        CigarParser(c, ref, b2).process(region, vd);
        if (!c.disableSV) filterSVStructures(vd, vd.maxReadLength);
        adjustMNP(vd, ref, c, region);
        auto reload = [&](int ms, int me) {
            Region rr; rr.chr = region.chr; rr.start = ms - 200; rr.end = me + 200; rr.gene = region.gene;
            CigarParser(c, ref, b1).process(rr, vd, /*reloadMode=*/true);
            CigarParser(c, ref, b2).process(rr, vd, /*reloadMode=*/true);
        };
        realigndel(vd, ref, c, region, vd.maxReadLength, {&b1, &b2});
        realignins(vd, ref, c, region, vd.maxReadLength);
        realignlgdel(vd, ref, c, region, vd.maxReadLength, reload, {&b1, &b2});
        realignlgins30(vd, ref, c, region, vd.maxReadLength, {&b1, &b2});
        realignlgins(vd, ref, c, region, vd.maxReadLength, reload, {&b1, &b2});
        if (!c.disableSV) findDEL(vd, ref, c, region, vd.maxReadLength, reload);
        if (!c.disableSV) findINV(vd, ref, c, region, vd.maxReadLength, reload, {&b1, &b2});
        if (!c.disableSV) findsv(vd, ref, c, region, vd.maxReadLength);
        if (!c.disableSV) findDELdisc(vd, ref, c, region, vd.maxReadLength);
        adjSNV(vd, ref);
        return vd;
    };

    auto processRegion = [&](const Region& region, Reference& ref, BamReader& bam, BamReader* bam2) {
        std::string buf;
        if (c.somatic) {
            // SomaticMode: run the pipeline on the tumor (BAM1) and normal (BAM2) separately, then
            // SomaticPostProcessModule compares the two per-position candidate sets. pos1 is built while
            // the reference window is loaded for BAM1; BAM2's pipeline reloads the same window (identical
            // ref bases), so pos1/pos2 both see consistent reference context.
            VariationData vd1 = runPipeline(ref, bam, region);
            auto pos1 = callVariantsSomatic(c, region, vd1, ref);
            VariationData vd2 = runPipeline(ref, *bam2, region);
            auto pos2 = callVariantsSomatic(c, region, vd2, ref);
            int maxRL = std::max(vd1.maxReadLength, vd2.maxReadLength);
            // combineAnalysis callback: re-run the merged pipeline over a widened window and return the
            // variant at (position, desc). pos1/pos2 are already built, so reloading `ref` here is safe.
            CombineFn combine = [&, maxRL](int vstart, int vend, int position,
                                           const std::string& desc, Variant& outv) -> bool {
                Region r2; r2.chr = region.chr; r2.start = vstart - maxRL; r2.end = vend + maxRL; r2.gene = region.gene;
                VariationData vdm = runMergedPipeline(ref, bam, *bam2, r2);
                auto posm = callVariantsSomatic(c, r2, vdm, ref);
                for (const auto& sp : posm) {
                    if (sp.position != position) continue;
                    for (const auto& vv : sp.variants)
                        if (vv.descriptionString == desc) { outv = vv; return true; }
                }
                return false;
            };
            appendSomaticRegion(buf, c, region, pos1, pos2, combine);
        } else {
            VariationData vd = runPipeline(ref, bam, region);
            auto variants = callVariants(c, region, vd, ref);
            for (const auto& v : variants) appendVariant(buf, c, region, v);
        }
        return buf;
    };

    const int nreg = (int)regions.size();
    int nthreads = std::max(1, std::min(c.threads, nreg));
    if (nthreads == 1) {
        Reference ref(c.ref);
        BamReader bam(c.bam);
        std::unique_ptr<BamReader> bam2;
        if (c.somatic) bam2.reset(new BamReader(c.bam2));
        for (const auto& region : regions) std::fputs(processRegion(region, ref, bam, bam2.get()).c_str(), stdout);
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
        BamReader bam(c.bam); // per-thread BAM file/index/header handle
        std::unique_ptr<BamReader> bam2;
        if (c.somatic) bam2.reset(new BamReader(c.bam2));
        int i;
        while ((i = nextWork.fetch_add(1)) < nreg) {
            std::string buf = processRegion(regions[i], ref, bam, bam2.get());
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

int main(int argc, char** argv) {
    // Convert any error (bad option value, unreadable/unindexed BAM, missing reference, malformed BED
    // or region, ...) into a single clear "vardictcpp: <what>" line and a non-zero exit, instead of an
    // uncaught-exception "terminate called" abort that tells the user nothing.
    try {
        return run(argc, argv);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "vardictcpp: %s\n", e.what());
        return 1;
    } catch (...) {
        std::fprintf(stderr, "vardictcpp: an unknown error occurred\n");
        return 1;
    }
}
