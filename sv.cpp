#include "sv.hpp"

using namespace std ;

SV::SV() { l = 0; }

SV::SV(const string type_, const string &chrom_, uint s_, const string &refall_, const string &altall_, const uint w_, const uint cov_, const int ngaps_, const int score_, bool imprecise_, uint l_) {
    type = type_ ;
    chrom = chrom_;
    s = s_;
    refall = refall_;
    altall = altall_;
    e = s + refall.size() - 1;
    w = w_;
    cov = cov_;
    ngaps = ngaps_;
    score = score_;
    imprecise = imprecise_;
    l = imprecise ? l_ : (int)altall.size() - (int)refall.size();
    idx = type + "_" + chrom + ":" + to_string(s) + "-" + to_string(e);
    idx += "_" + to_string(abs(l));
    gt = "./.";
}

void SV::genotype() {
    // gt = './.';
}

ostream &operator<<(ostream &os, const SV &sv) {
    os << sv.chrom << "\t"
       << sv.s << "\t"
       << sv.idx
       << "\t"
       << sv.refall << "\t"
       << sv.altall << "\t"
       << "."
       << "\t"
       << "PASS"
       << "\t"
       // INFO
       << "VARTYPE=SV;"
       << "SVTYPE=" << sv.type << ";"
       << "SVLEN=" << sv.l << ";"
       << "END=" << sv.e << ";"
       << "WEIGHT=" << sv.w << ";"
       << "COV=" << sv.cov << ";"
       << "AS=" << sv.score << ";"
       << "NV=" << sv.ngaps
       << (sv.imprecise ? ";IMPRECISE\t" : "\t")
       // -
       << "GT"
       << "\t"
       << sv.gt;
    return os;
}

void SVCluster::add_sv(SV sv, const int count) {
    svs[sv] += count ;
    s = min(s, int(sv.s)) ;
}
