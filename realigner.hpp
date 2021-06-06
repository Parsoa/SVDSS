#ifndef REALIGNER_HPP
#define REALIGNER_HPP 

#include <omp.h>
#include <vector>
#include <string>
#include <vector>
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

#include "bam.hpp"
#include "fastq.hpp"
#include "lprint.hpp"
#include "config.hpp"
#include "sfsutils.hpp"

using namespace std;

class Realigner {

public:
    void run();

private:

    int current_batch = 0 ;
    int current_input_batch = 0 ;
    int last_dumped_batch = 0 ;

    samFile *bam_file ;
    bam_hdr_t *bam_header ;
    std::vector<std::vector<std::vector<bam1_t*>>> bam_entries ;
    bool load_batch_bam(int threads, int batch_size, int p) ;
    void load_input_sfs_batch() ;
    ofstream out_file ;

    std::vector<std::unordered_map<std::string, std::vector<SFS>>> sfs_batches ;
    std::vector<std::string> process_batch(std::vector<bam1_t*> &bam_entries) ;
    std::vector<std::vector<std::vector<std::string>>> batches ;
    void output_batch(int) ;

    Configuration* config ;

    CIGAR rebuild_cigar(const std::string &ref_seq, const std::string &read_seq, const std::vector<std::pair<int, int>> &alpairs) ;
};

#endif
