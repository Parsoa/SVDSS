#ifndef INT_HPP
#define INT_HPP

#include <string>
#include <unordered_map>

#include "bed_utils.hpp"

class Intersector {

public:

    void load_tracks() ;

    std::unordered_map<std::string, std::vector<Track>> tracks ;

    Track* find(std::string chrom, int position, int len) ;

} ;

#endif
