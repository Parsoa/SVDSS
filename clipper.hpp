#ifndef CLIPPER_HPP
#define CLIPPER_HPP

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

struct Clip {
  std::string name;
  std::string chrom;
  uint p;
  uint l;
  bool starting;
  uint w;

  Clip() { w = 0; }

  Clip(std::string name_, std::string chrom_, const uint p_, uint l_,
       bool starting_, uint w_ = 0) {
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
  std::vector<Clip> clips;
  std::vector<Clip> remove_duplicates(const std::vector<Clip> &);
  std::vector<Clip> combine(const std::vector<Clip> &);
  std::vector<Clip> filter_lowcovered(const std::vector<Clip> &, const uint);
  std::vector<Clip> cluster(const std::vector<Clip> &, uint);
  std::vector<Clip>
  filter_tooclose_clips(const std::vector<Clip> &,
                        lib_interval_tree::interval_tree_t<int> &);

public:
  std::vector<std::vector<SV>> _p_svs;

  Clipper(const std::vector<Clip> &);
  void call(int threads, lib_interval_tree::interval_tree_t<int> &);
};

#endif
