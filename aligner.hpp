#ifndef ALIGNER_HPP
#define ALIGNER_HPP

#include <fstream>
#include <iostream>
#include <list>
#include <string>
#include <vector>
#include <zlib.h>

#include "config.hpp"
#include "bam.hpp"
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "kseq.h"
#include "sfsutils.hpp"

using namespace std;

class Aligner {
public:
  void run();

private:
  vector<pair<int, int>> get_aligned_pairs(bam1_t *aln);
  CIGAR rebuild_cigar(const string &chr, const string &qseq,
                      const vector<pair<int, int>> &alpairs);
};

#endif
