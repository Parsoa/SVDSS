#ifndef EXTENDER_HPP
#define EXTENDER_HPP

#include <ctime>
#include <fstream>
#include <omp.h>
#include <set>
#include <stdlib.h>
#include <string>
#include <vector>

#include <htslib/hts.h>
#include <htslib/sam.h>
// #include <interval_tree.hpp>
#include <parasail.h>
#include <parasail/matrices/nuc44.h>
#include <rapidfuzz/fuzz.hpp>
#include <spdlog/spdlog.h>

#include "bam.hpp"
#include "chromosomes.hpp"
#include "clipper.hpp"
#include "cluster.hpp"
#include "config.hpp"
#include "sfs.hpp"
#include "sv.hpp"

using namespace std;
// using namespace lib_interval_tree;

class Extender {

private:
  // book keeping:
  uint unplaced = 0; // SFS skipped since no first/last base can be placed from
                     // read alignment (should be rare)
  uint s_unplaced = 0;
  uint e_unplaced = 0;
  uint unknown = 0;
  uint unextended = 0;
  uint small_clusters =
      0; // number of cluster (before clustering) with low support
  uint small_extclusters =
      0; // number of extended clusters (after clustering) with low support

  Configuration *config;

  samFile *bam_file;
  hts_idx_t *bam_index;
  bam_hdr_t *bam_header;

  vector<Cluster> clusters;
  vector<ExtCluster> ext_clusters;
  unordered_map<string, vector<SFS>> *SFSs;
  vector<ExtSFS> extended_sfs;

  void extend_parallel();
  void extend_alignment(bam1_t *aln, int index);
  void extend_batch(vector<bam1_t *> bam_entries, int);
  bool load_batch_bam(int p);
  pair<int, int> get_unique_kmers(const vector<pair<int, int>> &alpairs,
                                  const uint k, const bool from_end,
                                  string chrom);

  void extract_sfs_sequences();
  void cluster_interval_tree();
  void cluster_no_interval_tree();

  void call();
  void filter_sv_chains();
  vector<pair<uint, char>> parse_cigar(string);
  vector<Cluster> cluster_by_length(const Cluster &cluster);

  // parallelize
  vector<vector<SV>> _p_svs;
  vector<vector<Clip>> _p_clips;
  vector<vector<ExtSFS>> _p_extended_sfs;
  vector<vector<Cluster>> _p_clusters;
  vector<vector<Consensus>> _p_alignments;
  vector<vector<vector<bam1_t *>>> bam_entries;
  vector<map<pair<int, int>, vector<ExtSFS>>> _p_sfs_clusters;

public:
  Extender(unordered_map<string, vector<SFS>> *);

  vector<SV> svs;
  vector<Clip> clips;
  vector<Consensus> alignments;

  void run();
};

#endif
