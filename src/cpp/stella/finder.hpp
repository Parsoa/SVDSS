#ifndef FND_HPP
#define FND_HPP

#include <mutex>
#include <vector>
#include <string>
#include <iterator>
#include <unordered_map>

#include "counter.hpp"
#include "bed_utils.hpp"

class Finder {

public:

    void run() ;

private:

    void load_sequences() ;
    void dump_sequences() ;

    std::mutex cout_lock ;
    std::unordered_map<std::string, Locus> sequences ;

} ;

#endif
