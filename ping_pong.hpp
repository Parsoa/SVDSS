#ifndef PNG_HPP
#define PNG_HPP

#include <assert.h>
#include <chrono>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <pthread.h>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>
#include <zlib.h>

#include <mrope.h>
#include <rld0.h>
#include <rle.h>

#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <htslib/hts.h>
#include <htslib/hts_endian.h>
#include <htslib/sam.h>

#include "assembler.hpp"
#include "config.hpp"
#include "fastq.hpp"
#include "sfs.hpp"

using namespace std;

#define fm6_comp(a) ((a) >= 1 && (a) <= 4 ? 5 - (a) : (a))

#define fm6_set_intv(e, c, ik)                                                 \
  ((ik).x[0] = (e)->cnt[(int)(c)],                                             \
   (ik).x[2] = (e)->cnt[(int)(c) + 1] - (e)->cnt[(int)(c)],                    \
   (ik).x[1] = (e)->cnt[fm6_comp(c)], (ik).info = 0)

/** From ropebwt2 ********/
static unsigned char seq_nt6_table[128] = {
    0, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 1,
    5, 2, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 1, 5, 2, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5};

typedef SFS sfs_type_t;
typedef map<string, vector<sfs_type_t>> batch_type_t;

static const vector<string> int2char({"$", "A", "C", "G", "T", "N"});

class PingPong {

public:
  int index();
  int search();

  int num_output_batches;

private:
  Configuration *config;

  int bam_mode;
  int current_batch = 0;
  int last_dumped_batch = 0;
  int reads_processed = 0;
  int non_x_reads = 0;

  gzFile fastq_file;
  kseq_t *fastq_iterator;
  samFile *bam_file;
  bam_hdr_t *bam_header;

  unordered_map<string, bool> smoothed_reads;
  unordered_map<string, bool> ignored_reads;
  void load_smoothed_read_ids();

  vector<vector<vector<int>>> read_seq_lengths;
  vector<vector<vector<int>>> read_seq_max_lengths;
  vector<vector<vector<uint8_t *>>> read_seqs;
  vector<vector<vector<string>>> read_names;
  vector<vector<vector<bam1_t *>>> bam_entries;
  vector<vector<vector<fastq_entry_t>>> fastq_entries;
  bool load_batch_bam(int p);
  bool load_batch_fastq(int threads, int batch_size, int p);
  batch_type_t process_batch(rld_t *index, int p, int i);
  void ping_pong_search(rld_t *index, uint8_t *seq, int l,
                        vector<sfs_type_t> &solutions, int);
  void output_batch(int);

  vector<vector<batch_type_t>> obatches;

  bool check_solution(rld_t *index, string S);
  bool backward_search(rld_t *index, const uint8_t *P, int p2);
  fastq_entry_t get_solution(fastq_entry_t fqe, int s, int l);
};

#endif
