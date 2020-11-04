#ifndef AGG_HPP
#define AGG_HPP

#include <string>
#include <unordered_map>

class Aggregator {

public:

    void run() ;

    void load_sequences() ;
    void dump_sequences() ;

    int num_batches ;

    std::unordered_map<std::string, int> sequences ;

} ;

#endif
