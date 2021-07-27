#ifndef CALLER_HPP
#define CALLER_HPP

#include <iostream>
#include <list>

#include "htslib/sam.h"
#include "htslib/hts.h"

#include "config.hpp"
#include "chromosomes.hpp"
#include "vcf.hpp"
#include "sv.hpp"
#include "insdeller.hpp"
#include "clipler.hpp"

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