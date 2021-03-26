#ifndef SNP_HPP
#define SNP_HPP

#include <vector>
#include <string>
#include <unordered_map>

#include <zlib.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include "htslib/hfile.h"
#include "htslib/hts_endian.h"

#include "vcf.hpp"
#include "config.hpp"

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
    bool load_batch_bam(int threads, int batch_size, int p) ;
    std::unordered_map<vcf_variant_t, int> process_batch_1(std::vector<bam1_t*> bam_entries) ;
    int pass(int index) ;
    void stats() ;

    std::unordered_map<vcf_variant_t, int> confidence_scores ;
    // ordered list of variants
    std::unordered_map<std::string, std::vector<vcf_variant_t>> vcf_variants ;
    // pass functions (minimize code redundancy)
    std::unordered_map<vcf_variant_t, int> (* pass_functions [2])(std::vector<bam1_t*>) ;
    std::pair<int, int> find_variants_in_read(int pos, int len, std::string chrom) ;

};

#endif
