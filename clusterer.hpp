#ifndef CLUSTERER_HPP
#define CLUSTERER_HPP

#include <ctime>
#include <fstream>
#include <omp.h>
#include <set>
#include <stdlib.h>
#include <string>
#include <vector>

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <spdlog/spdlog.h>

#include "bam.hpp"
#include "chromosomes.hpp"
#include "clipper.hpp"
#include "config.hpp"
#include "sfs.hpp"

using namespace std;

struct SubRead {
  string name;
  string seq;
  int htag; // 0: no tag; 1: hap1; 2: hap2

  SubRead(const string &_name, const string &_seq, const int _htag) {
    name = _name;
    seq = _seq;
    htag = _htag;
  }

  uint size() const { return seq.size(); }
};

struct Cluster {
  string chrom;
  int s;
  int e;
  int cov, cov0, cov1, cov2;
  vector<SFS> SFSs;
  vector<SubRead> subreads;

  Cluster() { cov = 0; };

  Cluster(const Cluster &c) {
    chrom = c.chrom;
    s = c.s;
    e = c.e;
    cov = c.cov;
    cov0 = c.cov0;
    cov1 = c.cov1;
    cov2 = c.cov2;
    SFSs = c.SFSs;
  };

  Cluster(const vector<SFS> &_SFSs) {
    SFSs = _SFSs;
    chrom = _SFSs[0].chrom;
  }

  Cluster(const string &_chrom, uint _s, uint _e, int _cov, int _cov0,
          int _cov1, int _cov2) {
    chrom = _chrom;
    s = _s;
    e = _e;
    cov = _cov;
    cov0 = _cov0;
    cov1 = _cov1;
    cov2 = _cov2;
  }

  void clear() {
    SFSs.clear();
  }

  void set_coordinates(int _s, int _e) {
    s = _s;
    e = _e;
  }

  void set_cov(vector<int> coverages) {
    cov0 = coverages[0];
    cov1 = coverages[1];
    cov2 = coverages[2];
    cov = cov0 + cov1 + cov2;
  }

  void add_subread(const string &name, const string &seq, int htag) {
    subreads.push_back(SubRead(name, seq, htag));
  }
  void add_subread(const SubRead &sr) { subreads.push_back(sr); }

  int get_len() const {
    uint l = 0;
    uint n = 0;
    for (const auto &sr : subreads) {
      ++n;
      l += sr.size();
    }
    return l / n;
  }

  SubRead get_subread(const int i) const { return subreads.at(i); }
  vector<string> get_names() const {
    vector<string> names;
    for (const SubRead &sr : subreads)
      names.push_back(sr.name);
    return names;
  }
  vector<string> get_seqs() const {
    vector<string> seqs;
    for (const SubRead &sr : subreads)
      seqs.push_back(sr.seq);
    return seqs;
  }
  string get_name(const int i) const { return subreads.at(i).name; }
  string get_seq(const int i) const { return subreads.at(i).seq; }
  SFS get_sfs(const int i) const { return SFSs.at(i); }

  uint size() const { return subreads.size(); }

  void print() const {
    cerr << chrom << ":" << s << "-" << e << " " << size() << "/" << cov
         << endl;
    for (const SubRead &sr : subreads)
      cerr << sr.name << " " << sr.size() << " " << sr.htag << " " << sr.seq
           << endl;
  }
};

class Clusterer {

private:
  Configuration *config;
  unordered_map<string, vector<SFS>> *SFSs;
  vector<SFS> extended_SFSs;
  samFile *bam_file;
  hts_idx_t *bam_index;
  bam_hdr_t *bam_header;

  // book keeping:
  uint unplaced = 0;   // SFS skipped since no first/last bases can be placed
  uint s_unplaced = 0; // SFS skipped since no first base can be placed
  uint e_unplaced = 0; // SFS skipped since no last base can be placed
  uint unknown = 0;    // SFS skipped due to (still) unknown reason
  uint unextended = 0; // SFS skipped since cannot be extended using kmers
  uint small_clusters =
      0; // number of clusters before extension  with low support
  uint small_clusters_2 =
      0; // number of clusters after extension with low support

  void align_and_extend();
  bool load_batch(int);
  void process_batch(int, int);
  void extend_alignment(bam1_t *, int);
  pair<int, int> get_unique_kmers(const vector<pair<int, int>> &alpairs,
                                  const uint k, const bool from_end,
                                  string chrom);
  void cluster_by_proximity();
  void fill_clusters();
  void store_clusters();

  // parallelize
  vector<vector<Clip>> _p_clips;
  vector<vector<SFS>> _p_extended_sfs;
  vector<vector<vector<bam1_t *>>> bam_entries;
  vector<map<pair<int, int>, vector<SFS>>> _p_sfs_clusters;

public:
  Clusterer(unordered_map<string, vector<SFS>> *);

  vector<Cluster> clusters;
  vector<Clip> clips;

  void run();
};

#endif
