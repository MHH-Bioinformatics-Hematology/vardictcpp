#include "var2vcf.hpp"
#include "config.hpp"
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>

namespace vardict {

// ---- small helpers ----------------------------------------------------------
static std::vector<std::string> splitTab(const std::string& s) {
    std::vector<std::string> f; size_t i = 0;
    while (true) {
        size_t j = s.find('\t', i);
        if (j == std::string::npos) { f.push_back(s.substr(i)); break; }
        f.push_back(s.substr(i, j - i)); i = j + 1;
    }
    return f;
}
static std::string get(const std::vector<std::string>& a, size_t i) { return i < a.size() ? a[i] : std::string(); }
static double toD(const std::string& s) { return s.empty() ? 0.0 : std::atof(s.c_str()); }
static long   toL(const std::string& s) { return s.empty() ? 0   : std::atol(s.c_str()); }

// Perl-style number stringification for the rare transformed odds ratio (1/x). Perl prints a float with
// up to 15 significant digits, trailing zeros stripped. %.15g matches that for these magnitudes.
static std::string perlNum(double v) {
    char b[64]; std::snprintf(b, sizeof(b), "%.15g", v); return b;
}
// format a frequency threshold as Perl interpolates the passed -f value (e.g. 0.01 -> "0.01")
static std::string freqStr(double f) { char b[32]; std::snprintf(b, sizeof(b), "%g", f); return b; }

// var2vcf_valid.pl reorder(): numeric chromosomes ascending by numeric value, then X, Y, M/MT, then the
// remaining (non-numeric) contigs in first-seen order.
static std::vector<std::string> reorder(const std::vector<std::string>& seen) {
    std::vector<std::pair<long, std::string>> num;
    std::vector<std::string> rest;
    bool hasX = false, haschrX = false, hasY = false, haschrY = false, hasMT = false, haschrM = false;
    for (const auto& c : seen) {
        bool hasDigit = false; for (char ch : c) if (isdigit((unsigned char)ch)) hasDigit = true;
        if (hasDigit && c.find('_') == std::string::npos) {
            std::string t; for (char ch : c) if (isdigit((unsigned char)ch)) t += ch;
            num.push_back({std::atol(t.c_str()), c});
        } else if (c == "X") haschrX = false, hasX = true;      // handled below via flags
        else if (c == "chrX") haschrX = true;
        else if (c == "Y") hasY = true;
        else if (c == "chrY") haschrY = true;
        else if (c == "MT") hasMT = true;
        else if (c == "chrM") haschrM = true;
        else rest.push_back(c);
    }
    std::stable_sort(num.begin(), num.end(), [](auto& a, auto& b) { return a.first < b.first; });
    std::vector<std::string> out;
    for (auto& p : num) out.push_back(p.second);
    if (hasX) out.push_back("X"); else if (haschrX) out.push_back("chrX");
    if (hasY) out.push_back("Y"); else if (haschrY) out.push_back("chrY");
    if (hasMT) out.push_back("MT"); else if (haschrM) out.push_back("chrM");
    for (auto& c : rest) out.push_back(c);
    return out;
}

static const char* HEADER_INFO_FILTER_FORMAT =
"##INFO=<ID=SAMPLE,Number=1,Type=String,Description=\"Sample name (with whitespace translated to underscores)\">\n"
"##INFO=<ID=TYPE,Number=1,Type=String,Description=\"Variant Type: SNV Insertion Deletion Complex\">\n"
"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n"
"##INFO=<ID=END,Number=1,Type=Integer,Description=\"Chr End Position\">\n"
"##INFO=<ID=VD,Number=1,Type=Integer,Description=\"Variant Depth\">\n"
"##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n"
"##INFO=<ID=BIAS,Number=1,Type=String,Description=\"Strand Bias Info\">\n"
"##INFO=<ID=REFBIAS,Number=1,Type=String,Description=\"Reference depth by strand\">\n"
"##INFO=<ID=VARBIAS,Number=1,Type=String,Description=\"Variant depth by strand\">\n"
"##INFO=<ID=PMEAN,Number=1,Type=Float,Description=\"The mean distance to the nearest 5 or 3 prime read end (whichever is closer) in all reads that support the variant call\">\n"
"##INFO=<ID=PSTD,Number=1,Type=Float,Description=\"Position STD in reads\">\n"
"##INFO=<ID=QUAL,Number=1,Type=Float,Description=\"Mean quality score in reads\">\n"
"##INFO=<ID=QSTD,Number=1,Type=Float,Description=\"Quality score STD in reads\">\n"
"##INFO=<ID=SBF,Number=1,Type=Float,Description=\"Strand Bias Fisher p-value\">\n"
"##INFO=<ID=ODDRATIO,Number=1,Type=Float,Description=\"Strand Bias Odds ratio\">\n"
"##INFO=<ID=MQ,Number=1,Type=Float,Description=\"Mean Mapping Quality\">\n"
"##INFO=<ID=SN,Number=1,Type=Float,Description=\"Signal to noise\">\n"
"##INFO=<ID=HIAF,Number=1,Type=Float,Description=\"Allele frequency using only high quality bases\">\n"
"##INFO=<ID=ADJAF,Number=1,Type=Float,Description=\"Adjusted AF for indels due to local realignment\">\n"
"##INFO=<ID=SHIFT3,Number=1,Type=Integer,Description=\"No. of bases to be shifted to 3 prime for deletions due to alternative alignment\">\n"
"##INFO=<ID=MSI,Number=1,Type=Float,Description=\"MicroSatellite. > 1 indicates MSI\">\n"
"##INFO=<ID=MSILEN,Number=1,Type=Float,Description=\"MicroSatellite unit length in bp\">\n"
"##INFO=<ID=NM,Number=1,Type=Float,Description=\"Mean mismatches in reads\">\n"
"##INFO=<ID=LSEQ,Number=1,Type=String,Description=\"5' flanking seq\">\n"
"##INFO=<ID=RSEQ,Number=1,Type=String,Description=\"3' flanking seq\">\n"
"##INFO=<ID=GDAMP,Number=1,Type=Integer,Description=\"No. of amplicons supporting variant\">\n"
"##INFO=<ID=TLAMP,Number=1,Type=Integer,Description=\"Total of amplicons covering variant\">\n"
"##INFO=<ID=NCAMP,Number=1,Type=Integer,Description=\"No. of amplicons don't work\">\n"
"##INFO=<ID=AMPFLAG,Number=1,Type=Integer,Description=\"Top variant in amplicons don't match\">\n"
"##INFO=<ID=HICNT,Number=1,Type=Integer,Description=\"High quality variant reads\">\n"
"##INFO=<ID=HICOV,Number=1,Type=Integer,Description=\"High quality total reads\">\n"
"##INFO=<ID=SPLITREAD,Number=1,Type=Integer,Description=\"No. of split reads supporting SV\">\n"
"##INFO=<ID=SPANPAIR,Number=1,Type=Integer,Description=\"No. of pairs supporting SV\">\n"
"##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"SV type: INV DUP DEL INS FUS\">\n"
"##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"The length of SV in bp\">\n"
"##INFO=<ID=DUPRATE,Number=1,Type=Float,Description=\"Duplication rate in fraction\">\n";

// FILTER + FORMAT block, thresholds interpolated (defaults from var2vcf_valid.pl; f uses cfg.freq).
static std::string headerFilters(double freq) {
    std::string f = freqStr(freq);
    return
"##FILTER=<ID=q22.5,Description=\"Mean Base Quality Below 22.5\">\n"
"##FILTER=<ID=Q10,Description=\"Mean Mapping Quality Below 10\">\n"
"##FILTER=<ID=p8,Description=\"Mean Position in Reads Less than 8\">\n"
"##FILTER=<ID=SN1.5,Description=\"Signal to Noise Less than 1.5\">\n"
"##FILTER=<ID=Bias,Description=\"Strand Bias\">\n"
"##FILTER=<ID=pSTD,Description=\"Position in Reads has STD of 0\">\n"
"##FILTER=<ID=d3,Description=\"Total Depth < 3\">\n"
"##FILTER=<ID=v2,Description=\"Var Depth < 2\">\n"
"##FILTER=<ID=f" + f + ",Description=\"Allele frequency < " + f + "\">\n"
"##FILTER=<ID=MSI12,Description=\"Variant in MSI region with 12 non-monomer MSI or 13 monomer MSI\">\n"
"##FILTER=<ID=NM5.25,Description=\"Mean mismatches in reads >= 5.25, thus likely false positive\">\n"
"##FILTER=<ID=InGap,Description=\"The variant is in the deletion gap, thus likely false positive\">\n"
"##FILTER=<ID=InIns,Description=\"The variant is adjacent to an insertion variant\">\n"
"##FILTER=<ID=Cluster0bp,Description=\"Two variants are within 0 bp\">\n"
"##FILTER=<ID=LongMSI,Description=\"The somatic variant is flanked by long A/T (>=14)\">\n"
"##FILTER=<ID=AMPBIAS,Description=\"Indicate the variant has amplicon bias.\">\n"
"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n"
"##FORMAT=<ID=VD,Number=1,Type=Integer,Description=\"Variant Depth\">\n"
"##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n"
"##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n"
"##FORMAT=<ID=RD,Number=2,Type=Integer,Description=\"Reference forward, reverse reads\">\n"
"##FORMAT=<ID=ALD,Number=2,Type=Integer,Description=\"Variant forward, reverse reads\">\n";
}

std::string var2vcfSingle(const Config& cfg, const std::string& tsv) {
    const int TotalDepth = 3, VarDepth = 2;
    const double Freq = cfg.freq, Pmean = 8, qmean = 22.5, Qmean = 10, GTFreq = 0.2, SN = 1.5;
    const int opt_I = 12; const double opt_m = 5.25; const int opt_c = 0, opt_P = 1, opt_T = 1;
    const bool noEnd = cfg.vcfNoEnd, passOnly = cfg.vcfPassOnly;

    // group rows by chr (first-seen order) -> pos -> rows
    std::vector<std::string> chrOrder;
    std::map<std::string, std::map<long, std::vector<std::vector<std::string>>>> hash;
    std::string sample;
    size_t i = 0;
    while (i < tsv.size()) {
        size_t nl = tsv.find('\n', i);
        std::string line = tsv.substr(i, nl == std::string::npos ? std::string::npos : nl - i);
        i = (nl == std::string::npos) ? tsv.size() : nl + 1;
        if (line.empty()) continue;
        if (line.find("R_HOME") != std::string::npos) continue;
        auto a = splitTab(line);
        if (a.size() < 8) continue;
        sample = a[0];
        const std::string& chr = a[2];
        if (hash.find(chr) == hash.end()) chrOrder.push_back(chr);
        hash[chr][toL(a[3])].push_back(std::move(a));
    }
    if (!cfg.sample.empty()) sample = cfg.sample;
    std::string sampleNW = sample; for (char& c : sampleNW) if (isspace((unsigned char)c)) c = '_';

    std::string out = "##fileformat=VCFv4.2\n##source=VarDict_v1.8.2\n";
    out += HEADER_INFO_FILTER_FORMAT;
    out += headerFilters(Freq);
    out += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample + "\n";
    if (hash.empty()) return out;

    auto chrs = reorder(chrOrder);
    for (const auto& chr : chrs) {
        auto cit = hash.find(chr); if (cit == hash.end()) continue;
        long pvs = 0;                     // previous PASS SNV start (Cluster filter)
        for (auto& pe : cit->second) {
            auto tmp = pe.second;
            std::stable_sort(tmp.begin(), tmp.end(), [](const std::vector<std::string>& x, const std::vector<std::string>& y) {
                return toD(get(y, 14)) < toD(get(x, 14));   // by AF descending
            });
            // default (no -A): only the top variant at the position
            const auto& a = tmp[0];
            std::string ref = get(a, 5); if (ref.empty()) continue;
            long start = toL(a[3]), end = toL(get(a, 4));
            std::string alt = get(a, 6), type = get(a, 35); if (type.empty()) type = "REF";
            long dp = toL(get(a, 7)); std::string vd = get(a, 8);
            long rfwd = toL(get(a, 9)), rrev = toL(get(a, 10)), vfwd = toL(get(a, 11)), vrev = toL(get(a, 12));
            std::string genotype = get(a, 13); double af = toD(get(a, 14));
            std::string bias = get(a, 15); std::string pmean = get(a, 16); std::string pstdS = get(a, 17);
            std::string qualS = get(a, 18); std::string qstd = get(a, 19); std::string sbf = get(a, 20);
            std::string oddratioS = get(a, 21); std::string mapqS = get(a, 22); std::string snS = get(a, 23);
            double hiaf = toD(get(a, 24)); std::string adjaf = get(a, 25); std::string shift3 = get(a, 26);
            std::string msiS = get(a, 27); std::string msilenS = get(a, 28); std::string nmS = get(a, 29);
            long hicnt = toL(get(a, 30)); long hicov = toL(get(a, 31));
            std::string lseq = get(a, 32), rseq = get(a, 33);
            std::string dupOrGamp = get(a, 36), svOrTamp = get(a, 37);
            bool isamp = a.size() >= 40;   // ampflag present only in amplicon mode
            double pstd = toD(pstdS), qual = toD(qualS), mapq = toD(mapqS), sn = toD(snS);
            double msi = toD(msiS), msilen = toD(msilenS), nm = toD(nmS), sbfv = toD(sbf);
            long rd = rfwd + rrev;

            // odds ratio: Inf -> 0; (0,1) -> 1/x
            std::string oddratio = oddratioS;
            if (oddratioS == "Inf") oddratio = "0";
            else { double o = toD(oddratioS); if (o < 1 && o > 0) oddratio = perlNum(1.0 / o); }
            double oddratioN = (oddratio == "Inf") ? 0 : toD(oddratio);

            std::vector<std::string> filters;
            if (dp < TotalDepth && !(hicnt * hiaf >= 0.5)) filters.push_back("d3");
            if (hicnt < VarDepth && !(hicnt * hiaf >= 0.5)) filters.push_back("v2");
            if (af < Freq) filters.push_back("f" + freqStr(Freq));
            if (toD(pmean) < Pmean) filters.push_back("p8");
            if (opt_P && pstd == 0 && !isamp && af < 0.35) filters.push_back("pSTD");
            if (qual < qmean) filters.push_back("q22.5");
            if (mapq < Qmean && af < 0.8) filters.push_back("Q10");
            if (sn < SN) filters.push_back("SN1.5");
            if (nm > opt_m) filters.push_back("NM5.25");
            long rl = (long)ref.size(), al = (long)alt.size();
            if ((msi > opt_I && msilen > 1 && af < 0.2 && std::labs(rl - al) == (long)msilen) ||
                (msi >= 13 && msilen == 1 && af <= 0.275 && std::labs(rl - al) == (long)msilen))
                filters.push_back("MSI12");
            if (hiaf < 0.25 && bias == "2;1" && sbfv < 0.01 && (oddratioN > 5 || oddratioN == 0) && end - start < 100)
                filters.push_back("Bias");
            if (std::labs(rl - al) == (long)msilen) {
                if (hiaf <= 0.275 && msi >= 13) filters.push_back("LongMSI");
                else if (hiaf <= 0.2 && msi >= 8 && msilen > 1) filters.push_back("LongMSI");
            }
            if (type == "SNV" && filters.empty() && start - pvs < opt_c) filters.push_back("Cluster0bp");
            if (isamp) {
                long gamp = toL(dupOrGamp), tamp = toL(svOrTamp), ncamp = toL(get(a, 38)); std::string ampflag = get(a, 39);
                if ((gamp < tamp - ncamp) || (!ampflag.empty() && ampflag != "0")) filters.push_back("AMPBIAS");
            }
            std::string filter = filters.empty() ? "PASS" : [&]{ std::string s; for (size_t k=0;k<filters.size();++k){ if(k) s+=";"; s+=filters[k]; } return s; }();
            if (passOnly && filter != "PASS") continue;

            std::string biasColon = bias; { size_t sc = biasColon.find(';'); if (sc != std::string::npos) biasColon[sc] = ':'; }
            long QUAL = (toL(vd) <= 1) ? 0 : (long)(std::log((double)toL(vd)) / std::log(2.0) * qual);
            std::string END = noEnd ? "" : (";END=" + std::to_string(end));

            // SV info (non-amplicon)
            std::string SVINFO; long splitreads = 0, spanpairs = 0;
            if (!isamp) {
                bool tampHasDash = svOrTamp.find('-') != std::string::npos;
                if (tampHasDash) {
                    auto parts = splitTab(std::string()); // placeholder; parse by '-'
                    std::vector<std::string> p; size_t s0 = 0;
                    for (size_t k = 0; k <= svOrTamp.size(); ++k) if (k == svOrTamp.size() || svOrTamp[k] == '-') { p.push_back(svOrTamp.substr(s0, k - s0)); s0 = k + 1; }
                    if (p.size() >= 2) { splitreads = toL(p[0]); spanpairs = toL(p[1]); }
                }
                if (alt.find('<') != std::string::npos) {
                    if (!(splitreads >= opt_T)) continue;   // skip SV without split support
                    long svlen = end - start; if (alt.find("INV") != std::string::npos) svlen++;
                    std::string sv = ";SVTYPE=" + alt + ";SVLEN=" + std::to_string(svlen);
                    std::string cleaned; for (char c : sv) if (c != '<' && c != '>') cleaned += c; SVINFO = cleaned;
                }
                if (!svOrTamp.empty()) SVINFO += ";SPLITREAD=" + std::to_string(splitreads) + ";SPANPAIR=" + std::to_string(spanpairs);
            }

            std::string gt;
            if (ref == alt) { alt = "."; gt = "0/0"; }
            else gt = (1 - af < GTFreq) ? "1/1" : (af >= 0.5 ? "1/0" : (af >= Freq ? "0/1" : "0/0"));
            std::string ad = (vd == "0") ? std::to_string(rd) : (std::to_string(rd) + "," + vd);
            std::string ampinfo, dupinfo;
            if (isamp) ampinfo = ";GDAMP=" + dupOrGamp + ";TLAMP=" + svOrTamp + ";NCAMP=" + get(a, 38) + ";AMPFLAG=" + get(a, 39);
            else if (!dupOrGamp.empty()) dupinfo = ";DUPRATE=" + dupOrGamp;

            std::string afS = get(a, 14);
            out += chr + "\t" + std::to_string(start) + "\t.\t" + ref + "\t" + alt + "\t" + std::to_string(QUAL)
                 + "\t" + filter + "\t"
                 + "SAMPLE=" + sampleNW + ";TYPE=" + type + ";DP=" + std::to_string(dp) + END + ";VD=" + vd
                 + ";AF=" + afS + ";BIAS=" + biasColon + ";REFBIAS=" + std::to_string(rfwd) + ":" + std::to_string(rrev)
                 + ";VARBIAS=" + std::to_string(vfwd) + ":" + std::to_string(vrev) + ";PMEAN=" + pmean + ";PSTD=" + pstdS
                 + ";QUAL=" + qualS + ";QSTD=" + qstd + ";SBF=" + sbf + ";ODDRATIO=" + oddratio + ";MQ=" + mapqS
                 + ";SN=" + snS + ";HIAF=" + get(a, 24) + ";ADJAF=" + adjaf + ";SHIFT3=" + shift3 + ";MSI=" + msiS
                 + ";MSILEN=" + msilenS + ";NM=" + nmS + ";HICNT=" + std::to_string(hicnt) + ";HICOV=" + std::to_string(hicov)
                 + ";LSEQ=" + lseq + ";RSEQ=" + rseq + ampinfo + dupinfo + SVINFO
                 + "\tGT:DP:VD:AD:AF:RD:ALD\t"
                 + gt + ":" + std::to_string(dp) + ":" + vd + ":" + ad + ":" + afS + ":" + std::to_string(rfwd) + ","
                 + std::to_string(rrev) + ":" + std::to_string(vfwd) + "," + std::to_string(vrev) + "\n";

            if (type == "SNV" && filter == "PASS") pvs = start;
        }
    }
    return out;
}

// header for the paired (somatic) VCF; thresholds are var2vcf_paired.pl defaults, f uses cfg.freq.
static std::string pairedHeader(double freq) {
    std::string f = freqStr(freq);
    return
"##fileformat=VCFv4.2\n##source=VarDict_v1.8.2\n"
"##INFO=<ID=SAMPLE,Number=1,Type=String,Description=\"Sample name (with whitespace translated to underscores)\">\n"
"##INFO=<ID=TYPE,Number=1,Type=String,Description=\"Variant Type: SNV Insertion Deletion Complex\">\n"
"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n"
"##INFO=<ID=END,Number=1,Type=Integer,Description=\"Chr End Position\">\n"
"##INFO=<ID=VD,Number=1,Type=Integer,Description=\"Variant Depth\">\n"
"##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n"
"##INFO=<ID=SHIFT3,Number=1,Type=Integer,Description=\"No. of bases to be shifted to 3 prime for deletions due to alternative alignment\">\n"
"##INFO=<ID=MSI,Number=1,Type=Float,Description=\"MicroSatellite. > 1 indicates MSI\">\n"
"##INFO=<ID=MSILEN,Number=1,Type=Float,Description=\"MSI unit repeat length in bp\">\n"
"##INFO=<ID=SSF,Number=1,Type=Float,Description=\"P-value\">\n"
"##INFO=<ID=SOR,Number=1,Type=Float,Description=\"Odds ratio\">\n"
"##INFO=<ID=LSEQ,Number=1,Type=String,Description=\"5' flanking seq\">\n"
"##INFO=<ID=RSEQ,Number=1,Type=String,Description=\"3' flanking seq\">\n"
"##INFO=<ID=STATUS,Number=1,Type=String,Description=\"Somatic or germline status\">\n"
"##INFO=<ID=P0.01Likely,Number=0,Type=Flag,Description=\"Likely candidate but p-value > 0.01/5**vd2 (means the evidence in tumor sample might be weak, e.g. small diff in AF)\">\n"
"##INFO=<ID=InDelLikely,Number=0,Type=Flag,Description=\"Likely indels more than 2bp are not considered somatic (weak evidence of presence in normal samples)\">\n"
"##FILTER=<ID=q22.5,Description=\"Mean Base Quality Below 22.5\">\n"
"##FILTER=<ID=Q0,Description=\"Mean Mapping Quality Below 0\">\n"
"##FILTER=<ID=p8,Description=\"Mean Position in Reads Less than 8\">\n"
"##FILTER=<ID=SN1.5,Description=\"Signal to Noise Less than 1.5\">\n"
"##FILTER=<ID=Bias,Description=\"Strand Bias\">\n"
"##FILTER=<ID=pSTD,Description=\"Position in Reads has STD of 0\">\n"
"##FILTER=<ID=MAF0.05,Description=\"Matched sample has AF > 0.05, thus not somatic\">\n"
"##FILTER=<ID=d5,Description=\"Total Depth < 5\">\n"
"##FILTER=<ID=v3,Description=\"Var Depth < 3\">\n"
"##FILTER=<ID=f" + f + ",Description=\"Allele frequency < " + f + "\">\n"
"##FILTER=<ID=P0.05,Description=\"Not significant with p-value > 0.05\">\n"
"##FILTER=<ID=DIFF0.2,Description=\"Non-somatic or LOH and allele frequency difference < 0.2\">\n"
"##FILTER=<ID=MSI12,Description=\"Variant in MSI region with 12 non-monomer MSI or 12 monomer MSI\">\n"
"##FILTER=<ID=NM5.25,Description=\"Mean mismatches in reads >= 5.25, thus likely false positive\">\n"
"##FILTER=<ID=InGap,Description=\"The somatic variant is in the deletion gap, thus likely false positive\">\n"
"##FILTER=<ID=InIns,Description=\"The somatic variant is adjacent to an insertion variant\">\n"
"##FILTER=<ID=Cluster0bp,Description=\"Two somatic variants are within 0 bp\">\n"
"##FILTER=<ID=LongAT,Description=\"The somatic variant is flanked by long A/T (>=14)\">\n"
"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n"
"##FORMAT=<ID=VD,Number=1,Type=Integer,Description=\"Variant Depth\">\n"
"##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n"
"##FORMAT=<ID=ALD,Number=2,Type=Integer,Description=\"Variant forward, reverse reads\">\n"
"##FORMAT=<ID=RD,Number=2,Type=Integer,Description=\"Reference forward, reverse reads\">\n"
"##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n"
"##FORMAT=<ID=ADJAF,Number=1,Type=Float,Description=\"Adjusted AF for indels due to local realignment\">\n"
"##FORMAT=<ID=BIAS,Number=2,Type=String,Description=\"Strand Bias Info\">\n"
"##FORMAT=<ID=PMEAN,Number=1,Type=Float,Description=\"The mean distance to the nearest 5 or 3 prime read end (whichever is closer) in all reads that support the variant call\">\n"
"##FORMAT=<ID=PSTD,Number=1,Type=Float,Description=\"Position STD in reads\">\n"
"##FORMAT=<ID=QUAL,Number=1,Type=Float,Description=\"Mean quality score in reads\">\n"
"##FORMAT=<ID=QSTD,Number=1,Type=Float,Description=\"Quality score STD in reads\">\n"
"##FORMAT=<ID=SBF,Number=1,Type=Float,Description=\"Strand Bias Fisher p-value\">\n"
"##FORMAT=<ID=ODDRATIO,Number=1,Type=Float,Description=\"Strand Bias Odds ratio\">\n"
"##FORMAT=<ID=MQ,Number=1,Type=Integer,Description=\"Mean Mapping Quality\">\n"
"##FORMAT=<ID=SN,Number=1,Type=Float,Description=\"Signal to noise\">\n"
"##FORMAT=<ID=HIAF,Number=1,Type=Float,Description=\"Allele frequency using only high quality bases\">\n"
"##FORMAT=<ID=NM,Number=1,Type=Float,Description=\"Mean mismatches in reads\">\n";
}

static std::string round2(double v) { char b[32]; std::snprintf(b, sizeof(b), "%.2f", v); return b; }
static std::string round0(double v) { char b[32]; std::snprintf(b, sizeof(b), "%.0f", v); return b; }

std::string var2vcfPaired(const Config& cfg, const std::string& tsv) {
    const int MinDepth = 5, VarDepth = 3;
    const double FREQ = cfg.freq, PMEAN = 8, QMEAN = 22.5, MQMEAN = 0, GTFREQ = 0.2, SN = 1.5;
    const int opt_I = 12; const double opt_m = 5.25; const int opt_c = 0;
    const bool passOnly = cfg.vcfPassOnly;   // opt_M (somatic-only) stays off, matching the wrapper

    std::vector<std::string> chrOrder;
    std::map<std::string, std::map<long, std::vector<std::vector<std::string>>>> hash;
    std::string sample = "tumor";
    size_t i = 0;
    while (i < tsv.size()) {
        size_t nl = tsv.find('\n', i);
        std::string line = tsv.substr(i, nl == std::string::npos ? std::string::npos : nl - i);
        i = (nl == std::string::npos) ? tsv.size() : nl + 1;
        if (line.empty() || line.find("R_HOME") != std::string::npos) continue;
        auto a = splitTab(line); if (a.size() < 8) continue;
        sample = a[0];
        if (hash.find(a[2]) == hash.end()) chrOrder.push_back(a[2]);
        hash[a[2]][toL(a[3])].push_back(std::move(a));
    }
    std::string samplem = sample + "-match";
    if (!cfg.sample.empty()) {
        sample = cfg.sample;
        samplem = cfg.sample2.empty() ? sample + "-match" : cfg.sample2;
    }
    std::string sampleNW = sample; for (char& c : sampleNW) if (isspace((unsigned char)c)) c = '_';

    std::string out = pairedHeader(FREQ);
    out += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample + "\t" + samplem + "\n";
    if (hash.empty()) return out;

    for (const auto& chr : reorder(chrOrder)) {
        auto cit = hash.find(chr); if (cit == hash.end()) continue;
        for (auto& pe : cit->second) {
            auto tmp = pe.second;
            std::stable_sort(tmp.begin(), tmp.end(), [](const std::vector<std::string>& x, const std::vector<std::string>& y) {
                return toD(get(y, 14)) < toD(get(x, 14));   // by tumor AF descending
            });
            const auto& a = tmp[0];
            std::string ref = get(a, 5); if (ref.empty()) continue;
            long start = toL(a[3]), end = toL(get(a, 4));
            std::string alt = get(a, 6);
            long dp1 = toL(get(a, 7)); double af1 = toD(get(a, 14));
            long vd1 = toL(get(a, 8)), rfwd1 = toL(get(a, 9)), rrev1 = toL(get(a, 10)), vfwd1 = toL(get(a, 11)), vrev1 = toL(get(a, 12));
            std::string bias1 = get(a, 15); double pmean1 = toD(get(a, 16)), pstd1 = toD(get(a, 17)), qual1 = toD(get(a, 18));
            double mapq1 = toD(get(a, 20)), sn1 = toD(get(a, 21)); double nm1 = toD(get(a, 24)); double sbf1 = toD(get(a, 25));
            std::string oddratio1S = get(a, 26);
            long dp2 = toL(get(a, 27)), vd2 = toL(get(a, 28)), rfwd2 = toL(get(a, 29)), rrev2 = toL(get(a, 30)), vfwd2 = toL(get(a, 31)), vrev2 = toL(get(a, 32));
            double af2 = toD(get(a, 34)); std::string bias2 = get(a, 35); double pmean2 = toD(get(a, 36)), pstd2 = toD(get(a, 37)), qual2 = toD(get(a, 38));
            double mapq2 = toD(get(a, 40)), sn2 = toD(get(a, 41)); double nm2 = toD(get(a, 44)); double sbf2 = toD(get(a, 45));
            std::string oddratio2S = get(a, 46);
            std::string shift3 = get(a, 47), msiS = get(a, 48), msilenS = get(a, 49), lseq = get(a, 50), rseq = get(a, 51);
            std::string status = get(a, 53), type = get(a, 54); if (type.empty()) type = "REF";
            std::string pvalue = get(a, 59); std::string oddratioS = get(a, 60);
            double msi = toD(msiS), msilen = toD(msilenS);
            long rd1 = rfwd1 + rrev1, rd2 = rfwd2 + rrev2;

            std::string oddratio = (oddratioS == "Inf") ? "0" : oddratioS;
            std::string oddratio1 = oddratio1S; { if (oddratio1S == "Inf") oddratio1 = "0"; else { double o = toD(oddratio1S); if (o < 1 && o > 0) oddratio1 = round2(1.0 / o); } }
            std::string oddratio2 = oddratio2S; { if (oddratio2S == "Inf") oddratio2 = "0"; else { double o = toD(oddratio2S); if (o < 1 && o > 0) oddratio2 = round2(1.0 / o); } }
            double oddratio1N = (oddratio1 == "Inf") ? 0 : toD(oddratio1), oddratio2N = (oddratio2 == "Inf") ? 0 : toD(oddratio2);

            std::vector<std::string> F, F2;
            bool strong = (status == "StrongSomatic");
            double pv = toD(pvalue);
            if (dp1 < MinDepth && !(strong && pv < 0.15 && af1 * vd1 >= 0.5)) F.push_back("d5");
            if (vd1 < VarDepth && !(strong && pv < 0.15 && af1 * vd1 >= 0.5)) F.push_back("v3");
            if (dp2 < MinDepth) F2.push_back("d5");
            if (vd2 < VarDepth) F2.push_back("v3");
            if (af1 < FREQ) F.push_back("f" + freqStr(FREQ));
            if (pmean1 < PMEAN) F.push_back("p8");
            if (pstd1 == 0 && vd1 < MinDepth) F.push_back("pSTD");
            if (qual1 < QMEAN) F.push_back("q22.5");
            if (mapq1 < MQMEAN) F.push_back("Q0");
            if (mapq1 < 10 && type == "SNV") F.push_back("Q0");
            if (sn1 < SN) F.push_back("SN1.5");
            if (nm1 >= opt_m) F.push_back("NM5.25");
            if (bias1 == "2;1" && sbf1 < 0.01 && (oddratio1N > 5 || oddratio1N == 0) && end - start < 100) F.push_back("Bias");
            if (af2 < FREQ) F2.push_back("f" + freqStr(FREQ));
            if (pmean2 < PMEAN) F2.push_back("p8");
            if (pstd2 == 0 && vd2 < MinDepth) F2.push_back("pSTD");
            if (qual2 < QMEAN) F2.push_back("q22.5");
            if (mapq2 < MQMEAN) F2.push_back("Q0");
            if (sn2 < SN) F2.push_back("SN1.5");
            if (nm2 >= opt_m) F2.push_back("NM5.25");
            bool hasBias = std::find(F.begin(), F.end(), "Bias") != F.end();
            if (!hasBias && bias2 == "2;1" && sbf2 < 0.01 && (oddratio2N > 5 || oddratio2N == 0) && end - start < 100) F.push_back("Bias");
            long rl = (long)ref.size(), al = (long)alt.size();
            if ((msi > opt_I && msilen > 1) || (msi > 12 && msilen == 1)) {
                if (!(strong && pv < 0.0005)) F.push_back("MSI12");
            }
            bool hasMsi = false; for (auto& s : F) if (s.rfind("MSI12", 0) == 0) hasMsi = true;
            if (std::labs(rl - al) == (long)msilen && !hasMsi) {
                if ((msi > opt_I && msilen > 1 && af1 < 0.35 && af2 < 0.35) || (msi > 12 && msilen == 1 && af1 < 0.35 && af2 < 0.35))
                    F.push_back("MSI12");
            }
            std::string filter = F.empty() ? "PASS" : [&]{ std::string s; for (size_t k=0;k<F.size();++k){ if(k) s+=";"; s+=F[k]; } return s; }();
            // germline rescue (opt_M off): a variant filtered only in tumor but clean in normal is PASS
            if (filter != "PASS" && F2.empty()) filter = "PASS";
            if (passOnly && filter != "PASS") continue;

            std::string gt = (1 - af1 < GTFREQ) ? "1/1" : (af1 >= 0.5 ? "1/0" : (af1 >= FREQ ? "0/1" : "0/0"));
            std::string gtm = (1 - af2 < GTFREQ) ? "1/1" : (af2 >= 0.5 ? "1/0" : (af2 >= FREQ ? "0/1" : "0/0"));
            std::string b1 = bias1; { size_t sc = b1.find(';'); if (sc != std::string::npos) b1[sc] = ','; } if (b1 == "0") b1 = "0,0";
            std::string b2 = bias2; { size_t sc = b2.find(';'); if (sc != std::string::npos) b2[sc] = ','; } if (b2 == "0") b2 = "0,0";
            std::string mq1 = round0(mapq1), mq2 = round0(mapq2);
            long QUAL = (vd1 > vd2) ? (long)(std::log((double)vd1) / std::log(2.0) * qual1)
                                    : (long)(std::log((double)vd2) / std::log(2.0) * qual2);

            std::string af1S = get(a, 14), af2S = get(a, 34), vd1S = get(a, 8), vd2S = get(a, 28);
            out += chr + "\t" + std::to_string(start) + "\t.\t" + ref + "\t" + alt + "\t" + std::to_string(QUAL) + "\t" + filter + "\t"
                 + "STATUS=" + status + ";SAMPLE=" + sampleNW + ";TYPE=" + type + ";DP=" + std::to_string(dp1) + ";VD=" + vd1S
                 + ";AF=" + af1S + ";SHIFT3=" + shift3 + ";MSI=" + msiS + ";MSILEN=" + msilenS + ";SSF=" + pvalue + ";SOR=" + oddratio
                 + ";LSEQ=" + lseq + ";RSEQ=" + rseq
                 + "\tGT:DP:VD:ALD:RD:AD:AF:BIAS:PMEAN:PSTD:QUAL:QSTD:SBF:ODDRATIO:MQ:SN:HIAF:ADJAF:NM\t"
                 + gt + ":" + std::to_string(dp1) + ":" + vd1S + ":" + std::to_string(vfwd1) + "," + std::to_string(vrev1)
                 + ":" + std::to_string(rfwd1) + "," + std::to_string(rrev1) + ":" + std::to_string(rd1) + "," + vd1S + ":" + af1S
                 + ":" + b1 + ":" + get(a, 16) + ":" + get(a, 17) + ":" + get(a, 18) + ":" + get(a, 19) + ":" + get(a, 25)
                 + ":" + oddratio1 + ":" + mq1 + ":" + get(a, 21) + ":" + get(a, 22) + ":" + get(a, 23) + ":" + get(a, 24) + "\t"
                 + gtm + ":" + std::to_string(dp2) + ":" + vd2S + ":" + std::to_string(vfwd2) + "," + std::to_string(vrev2)
                 + ":" + std::to_string(rfwd2) + "," + std::to_string(rrev2) + ":" + std::to_string(rd2) + "," + vd2S + ":" + af2S
                 + ":" + b2 + ":" + get(a, 36) + ":" + get(a, 37) + ":" + get(a, 38) + ":" + get(a, 39) + ":" + get(a, 45)
                 + ":" + oddratio2 + ":" + mq2 + ":" + get(a, 41) + ":" + get(a, 42) + ":" + get(a, 43) + ":" + get(a, 44) + "\n";
        }
    }
    return out;
}

} // namespace vardict
