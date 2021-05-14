#ifndef SNP_HPP
#define SNP_HPP

#include <vector>
#include <string>
#include <unordered_map>
#include <omp.h>
#include <ctime>
#include <chrono>
#include <thread>
#include <fstream>
#include <iomanip>
#include <assert.h>
#include <iostream>
#include <pthread.h>
#include <zlib.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <htslib/hts_endian.h>

#include "chromosomes.hpp"
#include "vcf.hpp"
#include "fastq.hpp"
#include "config.hpp"
#include "lprint.hpp"

using namespace std ;

#define BAM_CIGAR_MASK  0xf
#define BAM_CIGAR_TYPE  0x3C1A7

//TODO: why did I have to redefine this?
#define bam_set_seqi(s,i,b) ((s)[(i)>>1] = ((s)[(i)>>1] & (0xf0 >> ((~(i)&1)<<2))) | ((b)<<((~(i)&1)<<2)))

class SnpCorrector {

public:

    void run() ;

    // Loading BAM file
    samFile *bam_file ;
    bam_hdr_t* bam_header ;
    // output BAM file
    samFile* out_bam_file ;
    // <time < batch < reads > > >
    std::vector<std::vector<std::vector<bam1_t*>>> bam_entries ;

    // variant processing
    int correct_reads() ;
    bool load_batch_bam(int threads, int batch_size, int p) ;
    std::vector<fastq_entry_t> process_batch(std::vector<bam1_t*> bam_entries) ;
    fastq_entry_t correct_read(bam1_t* alignment, char* read_seq, std::string chrom) ;

    //std::unordered_map<vcf_variant_t, int> confidence_scores ;
    // ordered list of variants
    //std::unordered_map<std::string, std::vector<vcf_variant_t>> vcf_variants ;
    // pass functions (minimize code redundancy)
    //std::unordered_map<vcf_variant_t, int> (* pass_functions [2])(std::vector<bam1_t*>) ;
    //std::pair<int, int> find_variants_in_read(int pos, int len, std::string chrom) ;
    //void correct_snps(bam1_t* alignment, std::pair<int, int> limits, char* read_seq, std::string chrom) ;

};

#endif
