#ifndef SNP_HPP
#define SNP_HPP

#include <assert.h>
#include <chrono>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <pthread.h>
#include <set>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

#include <zlib.h>

#include <spdlog/spdlog.h>

#include "bam.hpp"
#include "chromosomes.hpp"
#include "config.hpp"
#include "fastq.hpp"

using namespace std;

class Smoother {

public:
  // Loading BAM file
  samFile *bam_file;
  bam_hdr_t *bam_header;
  hts_idx_t *bam_index;
  // output BAM file
  samFile *out_bam_file;
  vector<vector<vector<bam1_t *>>> bam_entries;

  vector<vector<vector<int>>> read_seq_max_lengths;
  vector<vector<vector<char *>>> read_seqs;
  vector<vector<vector<char *>>> new_read_seqs;
  vector<vector<vector<uint8_t *>>> new_read_quals;
  vector<vector<vector<int>>> new_read_seq_max_lengths;
  set<string> warned_chromosomes;
  void run();
  bool load_batch_bam(int p);
  void process_batch(vector<bam1_t *> bam_entries, int, int);
  void smooth_read(bam1_t *alignment, char *read_seq, int, int, int);

  Configuration *config;

private:
  double global_num_bases = 1;
  double global_num_mismatch;
  double global_num_indel;
  double expected_mismatch_rate = 0.002;
  int num_ignored_reads = 0;
  int reads_processed = 0;
};

#endif
