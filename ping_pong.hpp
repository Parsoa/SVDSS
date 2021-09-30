#ifndef PNG_HPP
#define PNG_HPP

#include <omp.h>
#include <ctime>
#include <vector>
#include <string>
#include <chrono>
#include <vector>
#include <thread>
#include <utility>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <assert.h>
#include <iostream>
#include <pthread.h>
#include <unordered_map>

#include "rle.h"
#include "rld0.h"
#include "mrope.h"

#include <zlib.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include "htslib/hfile.h"
#include "htslib/hts_endian.h"

#include "sfs.hpp"
#include "fastq.hpp"
#include "config.hpp"
#include "lprint.hpp"
#include "assembler.hpp"
#include "chromosomes.hpp"

using namespace std ;

#ifdef DEBUG_MODE
#  define DEBUG(x) x
#  define NEBUG(x)
#else
#  define DEBUG(x)
#  define NEBUG(x) x
#endif

#define fm6_comp(a) ((a) >= 1 && (a) <= 4? 5 - (a) : (a))

#define fm6_set_intv(e, c, ik) ((ik).x[0] = (e)->cnt[(int)(c)], (ik).x[2] = (e)->cnt[(int)(c)+1] - (e)->cnt[(int)(c)], (ik).x[1] = (e)->cnt[fm6_comp(c)], (ik).info = 0)

/** From ropebwt2 ********/
static unsigned char seq_nt6_table[128] = {
    0, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5
} ;

typedef SFS sfs_type_t ;
typedef std::map<std::string, std::vector<sfs_type_t>> batch_type_t; 

static const std::vector<std::string> int2char ({"$", "A", "C", "G", "T", "N"}) ;

class PingPong {

public:

    int index() ;
    int search() ;
    bool query(std::string) ;
    
    int num_output_batches ;

private:
    
    Configuration* config ;

    int mode ;
    int current_batch = 0 ;
    int last_dumped_batch = 0 ;

    gzFile fastq_file ;
    kseq_t* fastq_iterator ;
    samFile *bam_file ;
    bam_hdr_t *bam_header ;
    
    std::vector<std::vector<std::vector<int>>> read_seq_lengths ;
    std::vector<std::vector<std::vector<uint8_t*>>> read_seqs ;
    std::vector<std::vector<std::vector<bam1_t*>>> bam_entries ;
    std::vector<std::vector<std::vector<fastq_entry_t>>> fastq_entries ;
    bool load_batch_bam(int threads, int batch_size, int p) ;
    bool load_batch_fastq(int threads, int batch_size, int p) ;
    batch_type_t process_batch(rld_t* index, int p, int i, const bool isreversed) ;
    void ping_pong_search(rld_t *index, uint8_t* seq, int l, std::vector<sfs_type_t>& solutions, const bool isreversed) ;
    void output_batch(int) ;


    std::vector<std::vector<batch_type_t>> batches ;
    
    bool check_solution(rld_t* index, std::string S) ;
    bool backward_search(rld_t *index, const uint8_t *P, int p2) ;
    fastq_entry_t get_solution(fastq_entry_t fqe, int s, int l) ;

} ;

#endif
