#ifndef AGG_HPP
#define AGG_HPP

#include <string>
#include <vector>
#include <unordered_map>
#include <omp.h>
#include <iomanip>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include <string_view>

#include "chromosomes.hpp"
#include "lprint.hpp"

using namespace std ;

class Aggregator {

public:

    void run() ;

    void load_sequences() ;
    void dump_sequences() ;
    void find_high_abundance_sequences() ;

    int num_batches ;

    std::unordered_map<int, int> sequence_index ;
    std::unordered_map<std::string, int> sequences ;
    std::unordered_map<std::string, std::vector<std::string>> read_ids ;

} ;

#endif
