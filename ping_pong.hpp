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

extern "C" {
#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <htslib/hts.h>
#include <htslib/hts_endian.h>
#include <htslib/sam.h>
// fm-index must be included after htslib
#include <fm-index.h>
}
#include <spdlog/spdlog.h>

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

typedef map<string, vector<SFS>> batch_type_t;

static const vector<string> int2char({"$", "A", "C", "G", "T", "N"});

class PingPong {

public:
  int search();

private:
  Configuration *config;

  int bam_mode;
  int reads_processed = 0;

  gzFile fastq_file;
  kseq_t *fastq_iterator;
  samFile *bam_file;
  bam_hdr_t *bam_header;

  vector<vector<vector<uint>>> read_seq_lengths;
  vector<vector<vector<uint>>> read_seq_max_lengths;
  vector<vector<vector<uint8_t *>>> read_seqs;
  vector<vector<vector<string>>> read_names;
  vector<vector<vector<bam1_t *>>> bam_entries;
  vector<vector<vector<fastq_entry_t>>> fastq_entries;
  bool load_batch_bam(int);
  bool load_batch_fastq(int, int, int);
  batch_type_t process_batch(const rb3_fmi_t *, int, int);
  void ping_pong_search(const rb3_fmi_t *, const string &, uint8_t *, int,
                        vector<SFS> &, int);
  void output_batch(int);

  vector<vector<batch_type_t>> obatches;
};

#endif
