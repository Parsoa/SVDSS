#ifndef ALIGNER_HPP
#define ALIGNER_HPP

#include <fstream>
#include <iostream>
#include <list>
#include <string>
#include <vector>
#include <zlib.h>

#include "config.hpp"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "kseq.h"
#include "sfsutils.hpp"

using namespace std;

struct CIGAR {
  vector<int> ls;
  vector<char> ops;
  uint ed;

  CIGAR() { ed = 0; }

  void add(int l, char op, int e) {
    ed += e;
    if (ops.empty() || ops.back() != op) {
      ls.push_back(l);
      ops.push_back(op);
    } else {
      ls.back() += l;
    }
  }

  void add_front(int l) {
    ls[0] += l;
    ls.insert(ls.begin(), 1);
    ops.insert(ops.begin(), 'M');
  }

  void fixclips() {
    if (ops.front() != 'M')
      ops.front() = 'S';
    if (ops.back() != 'M')
      ops.back() = 'S';
  }

  string to_str() {
    string cigarstring;
    for (int i = 0; i < ls.size(); ++i) {
      cigarstring += to_string(ls[i]) + ops[i];
    }
    return cigarstring;
  }
};

class Aligner {
public:
  void run();

private:
  vector<pair<int, int>> get_aligned_pairs(bam1_t *aln);
  CIGAR rebuild_cigar(const string &chr, const string &qseq,
                      const vector<pair<int, int>> &alpairs);
};

#endif
