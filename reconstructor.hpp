#ifndef SNP_HPP
#define SNP_HPP

#include <vector>
#include <string>
#include <omp.h>
#include <ctime>
#include <chrono>
#include <thread>
#include <fstream>
#include <iomanip>
#include <assert.h>
#include <iostream>
#include <pthread.h>
#include <unordered_map>

#include <zlib.h>

#include "bam.hpp"
#include "vcf.hpp"
#include "fastq.hpp"
#include "config.hpp"
#include "lprint.hpp"
#include "chromosomes.hpp"

using namespace std ;

#define bam_set_seqi(s,i,b) ((s)[(i)>>1] = ((s)[(i)>>1] & (0xf0 >> ((~(i)&1)<<2))) | ((b)<<((~(i)&1)<<2)))

class Reconstructor {

public:

    // Loading BAM file
    samFile *bam_file ;
    bam_hdr_t* bam_header ;
    // output BAM file
    samFile* out_bam_file ;
    // <time < batch < reads > > >
    std::vector<std::vector<std::vector<bam1_t*>>> bam_entries ;

    void run() ;
    bool load_batch_bam(int threads, int batch_size, int p) ;
    void process_batch(std::vector<bam1_t*> bam_entries) ;
    void reconstruct_read(bam1_t* alignment, char* read_seq, std::string chrom) ;

    Configuration* config ;

private:

    double global_num_bases ;
    double global_num_mismatch ;
    double global_num_indel ;
    double expected_mismatch_rate = 0.002 ;
    int num_ignored_reads = 0 ;
};

#endif
