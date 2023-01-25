#include "sv.hpp"

using namespace std;

SV::SV() { l = 0; }

SV::SV(const string type_, const string &chrom_, uint s_, const string &refall_,
       const string &altall_, const uint w_, const uint cov_, const int ngaps_,
       const int score_, bool imprecise_, uint l_, string cigar_) {
  type = type_;
  chrom = chrom_;
  s = s_;
  refall = refall_;
  altall = altall_;
  e = s + refall.size() - 1;
  w = w_;
  l = l_;
  cov = cov_;
  ngaps = ngaps_;
  score = score_;
  imprecise = imprecise_;
  cigar = cigar_;
  idx = type + "_" + chrom + ":" + to_string(s) + "-" + to_string(e);
  idx += "_" + to_string(abs(l));
  gt = "./.";
  genotype();
}

void SV::add_reads(const vector<string> &names) {
  for (const string &name : names)
    reads += name + ",";
  reads.pop_back();
}

void SV::genotype() {
  if (imprecise) {
    gt = "./.";
  } else {
    float p = float(w) / float(cov);
    if (p >= 0.9)
      gt = "1/1";
    else
      gt = "0/1";
    // if (p <= 0.1)
    //     gt = "0/0";
  }
}

ostream &operator<<(ostream &os, const SV &sv) {
  os << sv.chrom << "\t" << sv.s << "\t" << sv.idx << "\t" << sv.refall << "\t"
     << sv.altall << "\t"
     << "."
     << "\t"
     << "PASS"
     << "\t"
     // INFO
     << "VARTYPE=SV;"
     << "SVTYPE=" << sv.type << ";"
     << "SVLEN=" << (sv.type == "DEL" ? -sv.l : sv.l) << ";"
     << "END=" << sv.e << ";"
     << "WEIGHT=" << sv.w << ";"
     << "COV=" << sv.cov << ";"
     << "AS=" << sv.score << ";"
     << "NV=" << sv.ngaps << ";"
     << "CIGAR=" << sv.cigar << ";"
     << "READS=" << sv.reads
     << (sv.imprecise ? ";IMPRECISE\t" : "\t")
     // -
     << "GT"
     << "\t" << sv.gt;
  return os;
}
