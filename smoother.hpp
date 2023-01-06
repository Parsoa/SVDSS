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
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

#include <zlib.h>

#include "bam.hpp"
#include "chromosomes.hpp"
#include "config.hpp"
#include "fastq.hpp"
#include "lprint.hpp"
#include "vcf.hpp"

using namespace std;

#define bam_set_seqi(s, i, b)                                                  \
  ((s)[(i) >> 1] =                                                             \
       ((s)[(i) >> 1] & (0xf0 >> ((~(i)&1) << 2))) | ((b) << ((~(i)&1) << 2)))

class Smoother {

public:
  // Loading BAM file
  samFile *bam_file;
  bam_hdr_t *bam_header;
  hts_idx_t *bam_index;
  // output BAM file
  samFile *out_bam_file;
  // <time < batch < reads > > >
  std::vector<std::vector<std::string>> smoothed_reads;
  std::vector<std::vector<std::string>> ignored_reads;
  std::vector<std::vector<std::vector<bam1_t *>>> bam_entries;

  std::vector<std::vector<std::vector<int>>> read_seq_lengths;
  std::vector<std::vector<std::vector<int>>> read_seq_max_lengths;
  std::vector<std::vector<std::vector<char *>>> read_seqs;
  std::vector<std::vector<std::vector<char *>>> new_read_seqs;
  std::vector<std::vector<std::vector<uint8_t *>>> new_read_quals;
  std::vector<std::vector<std::vector<int>>> new_read_seq_max_lengths;

  void run();
  bool load_batch_bam(int threads, int batch_size, int p);
  void process_batch(std::vector<bam1_t *> bam_entries, int, int);
  void smooth_read(bam1_t *alignment, char *read_seq, std::string chrom, int,
                   int, int);

  Configuration *config;

private:
  double global_num_bases = 1;
  double global_num_mismatch;
  double global_num_indel;
  double expected_mismatch_rate = 0.002;
  int num_ignored_reads = 0;
  int reads_processed = 0;

  void dump_smoothed_read_ids();
};

#endif
