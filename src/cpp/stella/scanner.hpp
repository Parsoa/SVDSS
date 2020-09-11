#ifndef SCN_HPP
#define SCN_HPP

#include <mutex>
#include <vector>
#include <string>
#include <iostream>
#include <iterator>

#include "counter.hpp"
#include "bed_utils.hpp"
#include "chromosomes.hpp"

typedef std::unordered_map<std::string, std::vector<Locus>> batch_type;

class Scanner {

public:

    void run() ;

private:
    
    void load_sequences() ;
    void scan_reference(int threads) ;
    void dump_sequences() ;
    batch_type scan_chromosome(std::string chrom, int, int) ;

    batch_type kmers ; // loci
    std::unordered_map<std::string, int> sequences ; // solutions

    std::mutex cout_mutex ;
} ;

#endif
