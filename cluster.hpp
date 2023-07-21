#ifndef CLUSTER_HPP
#define CLUSTER_HPP

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <abpoa.h>

#include "sfs.hpp"

using namespace std;

struct Cluster {
  string chrom;
  int s;
  int e;
  int cov;
  vector<string> names;
  vector<string> seqs;

  Cluster(const string &chrom_, uint s_, uint e_, uint cov_ = 0) {
    chrom = chrom_;
    s = s_;
    e = e_;
    cov = cov_;
  }

  void set_cov(uint cov_) { cov = cov_; }

  void add(const string &name, const string &seq) {
    names.push_back(name);
    seqs.push_back(seq);
  }

  int get_len() const {
    uint l = 0;
    uint n = 0;
    for (const string &seq : seqs) {
      ++n;
      l += seq.size();
    }
    return l / n;
  }

  string get_name(const int i) const { return names.at(i); }
  vector<string> get_names() const { return names; }
  string get_seq(const int i) const { return seqs.at(i); }

  uint size() const { return seqs.size(); }

  string poa();
};

struct ExtCluster {
  string chrom;
  vector<ExtSFS> seqs;

  ExtCluster(){};

  ExtCluster(const vector<ExtSFS> &_seqs) {
    seqs = _seqs;
    chrom = seqs[0].chrom;
  }
};

#endif
