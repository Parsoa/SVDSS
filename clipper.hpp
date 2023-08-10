#ifndef CLIPPER_HPP
#define CLIPPER_HPP

#include <algorithm>
#include <map>
#include <omp.h>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#include "htslib/hts.h"
#include "htslib/sam.h"
#include "interval_tree.hpp"

#include "chromosomes.hpp"
#include "sv.hpp"

using namespace std;
using namespace lib_interval_tree;

struct Clip {
  string name;
  string chrom;
  uint p;
  uint l;
  bool starting;
  uint w;

  Clip() { w = 0; }

  Clip(string name_, string chrom_, const uint p_, uint l_, bool starting_,
       uint w_ = 0) {
    name = name_;
    chrom = chrom_;
    p = p_;
    l = l_;
    starting = starting_;
    w = w_;
  }

  bool operator<(const Clip &c) const { return p < c.p; }
};

class Clipper {

private:
  vector<Clip> clips;
  vector<Clip> remove_duplicates(const vector<Clip> &);
  vector<Clip> combine(const vector<Clip> &);
  vector<Clip> filter_lowcovered(const vector<Clip> &, const uint);
  vector<Clip> cluster(const vector<Clip> &, uint);
  vector<Clip> filter_tooclose_clips(const vector<Clip> &,
                                     lib_interval_tree::interval_tree_t<int>
                                     &);

public:
  vector<vector<SV>> _p_svs;

  Clipper(const vector<Clip> &);
  void call(int threads, lib_interval_tree::interval_tree_t<int> &);
};

#endif
