#ifndef CALLER_HPP
#define CALLER_HPP

#include <iostream>

#include "htslib/sam.h"
#include "htslib/hts.h"

#include "sv.hpp"
#include "vcf.hpp"
#include "config.hpp"
#include "clipper.hpp"
#include "insdeller.hpp"
#include "chromosomes.hpp"

using namespace std;

class Caller {

public:
    void run();

private:
    Configuration* config;

    ofstream ovcf;
    ofstream osam;

    samFile *sfs_bam;
    bam_hdr_t *sfs_bamhdr;
    hts_idx_t *sfs_bamindex;

    samFile *read_bam;
    bam_hdr_t *read_bamhdr;
    hts_idx_t *read_bamindex;
};

#endif
