#ifndef INT_HPP
#define INT_HPP

#include <mutex>
#include <vector>
#include <string>
#include <iterator>
#include <unordered_map>

#include "counter.hpp"
#include "bed_utils.hpp"

class Intersector {

public:

    void run() ;

private:

    std::unordered_map<std::string, Locus> load_sequences(std::string) ;
    void dump_sequences() ;

    std::unordered_map<std::string, Locus> shared_sequences ;

} ;

#endif
