#ifndef SFS_HPP
#define SFS_HPP

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

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
  int s;
  int l;
  int htag;

  SFS() {
    s = 0;
    l = 0;
    htag = 0;
  }

  SFS(int s_, int l_, int htag_) {
    s = s_;
    l = l_;
    htag = htag_;
  }

  // void reverse(int p) { s = p - s - l; }
};

bool operator<(const SFS &, const SFS &);

struct ExtSFS {
  string chrom;
  string qname;
  // reference positions
  int s;
  int e;
  // query positions
  int qs;
  int qe;
  // other
  int htag;

  ExtSFS(const string &_chrom, const string &_qname, int _s, int _e, int _qs,
         int _qe, int _htag) {
    chrom = _chrom;
    qname = _qname;
    s = _s;
    e = _e;
    qs = _qs;
    qe = _qe;
    htag = _htag;
  }

  bool operator<(const ExtSFS &c) const {
    if (chrom == c.chrom) {
      return s < c.s;
    } else {
      return chrom < c.chrom;
    }
  }

  bool operator==(const ExtSFS &c) const {
    return chrom == c.chrom and s == c.s and e == c.e;
  }
};

class Consensus {
public:
  string seq;
  string chrom;
  string cigar;
  int s;
  int e;

  Consensus(const string _seq, const string _cigar, const string _chrom, int _s,
            int _e) {
    seq = _seq;
    cigar = _cigar;
    chrom = _chrom;
    s = _s;
    e = _e;
  }
};

map<string, vector<SFS>> parse_sfsfile(const string &, int);

#endif
