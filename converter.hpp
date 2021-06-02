#ifndef CONVERTER_HPP
#define CONVERTER_HPP

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
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include "htslib/hfile.h"
#include "htslib/hts_endian.h"

#include "fastq.hpp"
#include "lprint.hpp"
#include "config.hpp"
#include "sfsutils.hpp"

using namespace std;

class Converter {

public:
    void run();

private:

    int current_batch = 0 ;
    int current_input_batch = 0 ;

    gzFile fastq_file ;
    kseq_t* fastq_iterator ;
    samFile *bam_file ;
    bam_hdr_t *bam_header ;
    std::vector<std::vector<std::vector<bam1_t*>>> bam_entries ;
    std::vector<std::vector<std::vector<fastq_entry_t>>> fastq_entries ;
    bool load_batch_bam(int threads, int batch_size, int p) ;
    bool load_batch_fastq(int threads, int batch_size, int p) ;
    void load_input_sfs_batch() ;

    std::vector<std::unordered_map<std::string, std::vector<SFS>>> sfs_batches ;
    std::vector<fastq_entry_t> process_batch(std::vector<fastq_entry_t> &fastq_entries) ;
    std::vector<std::vector<std::vector<fastq_entry_t>>> batches ;
    void output_batch(int) ;

    Configuration* config ;
};

#endif
