#ifndef SFS_HPP
#define SFS_HPP

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include <spdlog/spdlog.h>

using namespace std;

static const char RCN[128] = {
    0,   0,   0, 0,   0,   0,   0,   0,   0,   0,   //  0
    0,   0,   0, 0,   0,   0,   0,   0,   0,   0,   // 10
    0,   0,   0, 0,   0,   0,   0,   0,   0,   0,   // 20
    0,   0,   0, 0,   0,   0,   0,   0,   0,   0,   // 30
    0,   0,   0, 0,   0,   0,   0,   0,   0,   0,   // 40
    0,   0,   0, 0,   0,   0,   0,   0,   0,   0,   // 50
    0,   0,   0, 0,   0,   'T', 0,   'G', 0,   0,   // 60
    0,   'C', 0, 0,   0,   0,   0,   0,   'N', 0,   // 70
    0,   0,   0, 0,   'A', 0,   0,   0,   0,   0,   // 80
    0,   0,   0, 0,   0,   0,   0,   'T', 0,   'G', // 90
    0,   0,   0, 'G', 0,   0,   0,   0,   0,   0,   // 100
    'N', 0,   0, 0,   0,   0,   'A', 0,   0,   0,   // 110
    0,   0,   0, 0,   0,   0,   0,   0              // 120
};

struct SFS {
  string chrom;
  string qname;
  // reference and query positions, length
  int rs;
  int re;
  int qs;
  int qe;
  int l;
  // other
  int htag; // 0: no tag; 1: hap1; 2: hap2

  SFS(const string &_qname, int _qs, int _l, int _htag) {
    chrom = "";
    qname = _qname;
    qs = _qs;
    qe = _qs + _l;
    l = _l;
    htag = _htag;
  }

  SFS(const string &_chrom, const string &_qname, int _rs, int _re, int _qs,
      int _qe, int _htag) {
    chrom = _chrom;
    qname = _qname;
    rs = _rs;
    re = _re;
    qs = _qs;
    qe = _qe;
    l = qe - qs + 1;
    htag = _htag;
  }

  // void reverse(int p) { s = p - s - l; }

  bool operator<(const SFS &c) const {
    // FIXME: bad trick to manage both cases we can have (in assembler, we have
    // to check for query positions and we are sure we do not have a chromosome
    // there. On caller, where we have the chrom, for reference positions)
    if (chrom == "" || c.chrom == "")
      return qs < c.qs;
    else
      return chrom == c.chrom ? rs < c.rs : chrom < c.chrom;
  }

  bool operator==(const SFS &c) const {
    return chrom == c.chrom and rs == c.rs and re == c.re;
  }
};

unordered_map<string, vector<SFS>> parse_sfsfile(const string &);

#endif
