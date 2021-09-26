#include "sv.hpp"

SV::SV() { l = 0; }

SV::SV(const string type_, const string &chrom_, uint s_, const string &refall_, const string &altall_, const uint w_, const uint cov_, int l_, int nv_, int score_, bool imprecise_)
{
    type = type_;
    chrom = chrom_;
    s = s_;
    refall = refall_;
    altall = altall_;
    e = s + refall.size() - 1;
    l = l_;
    w = w_;
    cov = cov_;
    score = score_;
    nv = nv_;
    imprecise = imprecise_;
    genotype();
}

void SV::genotype()
{
    if (imprecise)
        gt = "./.";
    else
    {
        float p = (float)w / (float)cov;
        if (p <= 0.1)
            gt = "0/0";
        else if (0.1 < p && p < 0.9)
            gt = "0/1";
        else
            gt = "1/1";
    }
}

ostream &operator<<(ostream &os, const SV &sv)
{
    os << sv.chrom << "\t"
       << sv.s << "\t"
       << "."
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
       << "NV=" << sv.nv << ";"
       << "ASCORE=" << sv.score
       << (sv.imprecise ? ";IMPRECISE" : "")
       << "\t"
       // -
       << "GT"
       << "\t"
       << sv.gt;
    return os;
}