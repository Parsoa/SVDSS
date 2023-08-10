#ifndef CALLER_HPP
#define CALLER_HPP

#include <ctime>
#include <dirent.h>
#include <iostream>
#include <map>

#include <abpoa.h>
#include <parasail.h>
#include <parasail/matrices/nuc44.h>
#include <rapidfuzz/fuzz.hpp>
#include <spdlog/spdlog.h>

#include "bam.hpp"
#include "chromosomes.hpp"
#include "clipper.hpp"
#include "clusterer.hpp"
#include "genotyper.hpp"
#include "config.hpp"
#include "sfs.hpp"
#include "sv.hpp"

using namespace std;

// AaCcGgTtNn ==> 0,1,2,3,4
static unsigned char _char26_table[256] = {
    0, 1,         2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4,         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4 /*'-'*/, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0,
    4, 1,         4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4,
    4, 4,         4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4,
    4, 4,         4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4,         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4,         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4,         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4,         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4,         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4,         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};

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

  friend ostream &operator<<(ostream &os, const Consensus &c) {
    os << c.chrom << ":" << c.s + 1 << "-" << c.e + 1 << "\t"
       << "0"
       << "\t" << c.chrom << "\t" << c.s + 1 << "\t"
       << "60"
       << "\t" << c.cigar << "\t"
       << "*"
       << "\t"
       << "0"
       << "\t"
       << "0"
       << "\t" << c.seq << "\t"
       << "*";
    return os;
  }
};

class Caller {

public:
  void run();

private:
  Configuration *config;
  unordered_map<string, vector<SFS>> SFSs;

  vector<SV> svs;
  vector<Consensus> alignments;

  void pcall(const vector<Cluster> &);
  void clean_dups();
  void filter_sv_chains();
  void write_vcf();
  void write_sam();

  vector<Cluster> split_cluster_by_len(const Cluster &);
  vector<Cluster> split_cluster(const Cluster &);
  string run_poa(const vector<string> &);

  // parallelize
  vector<vector<SV>> _p_svs;
  vector<vector<Consensus>> _p_alignments;

  void print_vcf_header();
};

#endif
