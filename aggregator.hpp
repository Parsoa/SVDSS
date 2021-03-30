#ifndef AGG_HPP
#define AGG_HPP

#include <string>
#include <vector>
#include <unordered_map>

class Aggregator {

public:

    void run() ;

    void load_sequences() ;
    void dump_sequences() ;
    void find_high_abundance_sequences() ;

    int num_batches ;

    std::unordered_map<int, int> sequence_index ;
    std::unordered_map<std::string, int> sequences ;
    std::unordered_map<std::string, std::unordered_map<std::string, int>> read_ids ;

} ;

#endif
