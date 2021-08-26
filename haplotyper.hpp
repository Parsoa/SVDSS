#ifndef HAPPER_HPP
#define HAPPER_HPP

#include <vector>
#include <string>
#include <assert.h>
#include <iostream>
#include <unordered_map>

#include "sv.hpp"
#include "bam.hpp"
#include "vcf.hpp"
#include "config.hpp"
#include "cluster.hpp"
#include "chromosomes.hpp"

using namespace std;

class Haplotyper {

public:

    Haplotyper() ;

    std::vector<SV> haplotype(Cluster cluster);

private:
    
    Configuration* config ;
    
    std::vector<SV> assemble_reads(Cluster cluster) ;

};

#endif
