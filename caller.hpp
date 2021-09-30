#ifndef CALLER_HPP
#define CALLER_HPP

#include <iostream>

#include "htslib/sam.h"
#include "htslib/hts.h"

#include "sv.hpp"
#include "vcf.hpp"
#include "config.hpp"
#include "clipper.hpp"
#include "extender.hpp"
#include "chromosomes.hpp"

using namespace std;

class Caller {

public:
    void run();

private:
    Configuration* config;
    std::unordered_map<std::string, std::vector<SFS>> SFSs ;

    ofstream ovcf;
    ofstream osam;
    void load_input_sfs() ;
    void print_vcf_header() ;
};

#endif
