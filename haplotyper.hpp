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
    std::vector<SV> map_to_chm13(Cluster cluster, std::string assembly) ;

private:
    
    Configuration* config ;
    std::vector<SV> call_svs(Cluster cluster, std::string bam_path, std::string ref_seq, std::string query, int offset) ;
    int make_working_dir(Cluster cluster) ;
    
    std::vector<SV> assemble_reads(Cluster cluster) ;

};

#endif
