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
  std::string type;
  std::string chrom;
  std::string idx;
  int s;
  int e;
  std::string refall;
  std::string altall;
  uint w;
  uint cov;
  int l = 0;
  int ngaps;
  int score;
  std::string gt;
  bool imprecise;
  string cigar;
  string reads;

  SV();
  SV(const std::string type_, const std::string &chrom_, uint s_,
     const std::string &refall_, const std::string &altall_, const uint w_,
     const uint cov_, const int ngaps_, const int score_,
     bool imprecise_ = false, uint l_ = 0, string cigar_ = ".");
  void add_reads(const vector<string> &reads_);
  void genotype();

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

  friend std::ostream &operator<<(std::ostream &os, const SV &sv);
};

#endif
