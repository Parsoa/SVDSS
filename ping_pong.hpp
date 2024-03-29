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

static inline int kputsn(const char *p, int l, kstring_t *s) {
  if (s->l + l + 1 >= s->m) {
    char *tmp;
    s->m = s->l + l + 2;
    kroundup32(s->m);
    if ((tmp = (char *)realloc(s->s, s->m)))
      s->s = tmp;
    else
      return EOF;
  }
  memcpy(s->s + s->l, p, l);
  s->l += l;
  s->s[s->l] = 0;
  return l;
}

typedef map<string, vector<SFS>> batch_type_t;

static const vector<string> int2char({"$", "A", "C", "G", "T", "N"});

class PingPong {

public:
  int index();
  int search();

private:
  Configuration *config;

  int bam_mode;
  int reads_processed = 0;

  gzFile fastq_file;
  kseq_t *fastq_iterator;
  samFile *bam_file;
  bam_hdr_t *bam_header;

  vector<vector<vector<int>>> read_seq_lengths;
  vector<vector<vector<int>>> read_seq_max_lengths;
  vector<vector<vector<uint8_t *>>> read_seqs;
  vector<vector<vector<string>>> read_names;
  vector<vector<vector<bam1_t *>>> bam_entries;
  vector<vector<vector<fastq_entry_t>>> fastq_entries;
  bool load_batch_bam(int);
  bool load_batch_fastq(int, int, int);
  batch_type_t process_batch(rld_t *, int, int);
  void ping_pong_search(rld_t *, const string &, uint8_t *, int, vector<SFS> &,
                        int);
  void output_batch(int);

  vector<vector<batch_type_t>> obatches;
};

#endif
