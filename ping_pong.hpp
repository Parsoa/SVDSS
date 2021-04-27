#ifndef PNG_HPP
#define PNG_HPP

#include <mutex>
#include <vector>
#include <string>
#include <iterator>
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
#include "fastq.hpp"
#include "config.hpp"

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

static const std::vector<std::string> int2char ({"$", "A", "C", "G", "T", "N"}) ;

struct OutputBatchArgs {
    int batch ;
} ;

class PingPong {

public:

    int index() ;
    bool query(std::string) ;
    int search() ;

    int num_output_batches ;

private:

    int current_batch = 0 ;

    gzFile fastq_file ;
    kseq_t* fastq_iterator ;
    std::vector<std::vector<std::vector<fastq_entry_t>>> fastq_entries ;
    bool load_batch_fastq(int threads, int batch_size, int p) ;
    std::vector<fastq_entry_t> process_batch_fastq(rld_t* index, std::vector<fastq_entry_t> fastq_entries) ;
 
    samFile *bam_file ;
    bam_hdr_t *bam_header ;
    std::vector<std::vector<std::vector<bam1_t*>>> bam_entries ;
    bool load_batch_bam(int threads, int batch_size, int p) ;
    std::vector<fastq_entry_t> process_batch_bam(rld_t* index, std::vector<bam1_t*> bam_entries) ;

    fastq_entry_t get_solution(fastq_entry_t fqe, int s, int l) ;

    std::vector<std::unordered_map<fastq_entry_t, int>> search_solutions ;
    std::vector<std::unordered_map<fastq_entry_t, std::vector<std::string>>> read_ids ;
    //void ping_pong_search(rld_t *index, const char* seq, const char* qual, std::vector<fastq_entry_t>& solutions) ;
    void ping_pong_search(rld_t *index, fastq_entry_t fqe, std::vector<fastq_entry_t>& solutions) ;
    void output_batch(void* args) ;
    
    bool backward_search(rld_t *index, const uint8_t *P, int p2) ;
    bool check_solution(rld_t* index, std::string S) ;

    Configuration* config ;
} ;

#endif
