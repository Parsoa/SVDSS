#ifndef SV_HPP
#define SV_HPP

#include <iostream>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

class SV {
public:
  string type;
  string chrom;
  string idx;
  int s;
  int e;
  string refall;
  string altall;
  uint w;
  int cov;
  int cov0;
  int cov1;
  int cov2;
  int l = 0;
  int ngaps;
  int score;
  string gt;
  bool imprecise;
  string cigar;
  string reads;

  SV();
  SV(const string type_, const string &chrom_, uint s_, const string &refall_,
     const string &altall_, const uint w_, const uint cov_, const int ngaps_,
     const int score_, bool imprecise_ = false, uint l_ = 0,
     string cigar_ = ".");
  void add_reads(const vector<string> &reads_);
  void genotype();

  void set_cov(int, int, int, int);

  bool operator<(const SV &c) const {
    if (chrom < c.chrom) {
      return true;
    } else if (chrom > c.chrom) {
      return false;
    } else {
      return s < c.s;
    }
  }

  bool operator==(const SV &c) const {
    return chrom == c.chrom and s == c.s and e == c.e and type == c.type;
  }

  friend ostream &operator<<(ostream &os, const SV &sv);
};

#endif
