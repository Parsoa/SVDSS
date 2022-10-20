#ifndef CLUSTER_HPP
#define CLUSTER_HPP

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "abpoa.h"

#include "sfs.hpp"

struct Cluster {
  std::string chrom;
  int s;
  int e;
  int cov;
  std::vector<std::string> seqs;

  Cluster(const std::string &chrom_, uint s_, uint e_, uint cov_ = 0) {
    chrom = chrom_;
    s = s_;
    e = e_;
    cov = cov_;
  }

  void set_cov(uint cov_) { cov = cov_; }

  void add(const std::string &seq) { seqs.push_back(seq); }

  int get_len() const {
    uint l = 0;
    uint n = 0;
    for (const std::string &seq : seqs) {
      ++n;
      l += seq.size();
    }
    return l / n;
  }

  std::vector<std::string> get_seqs() const { return seqs; }

  uint size() const { return seqs.size(); }

  void dump(std::ofstream &o) const;

  std::string poa();
};

struct ExtCluster {
  std::string chrom;
  std::vector<ExtSFS> seqs;

  ExtCluster(){};

  ExtCluster(const std::vector<ExtSFS> &_seqs) {
    seqs = _seqs;
    chrom = seqs[0].chrom;
  }
};

#endif
