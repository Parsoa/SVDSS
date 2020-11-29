#ifndef HPS_HPP
#define HPS_HPP

#include <string>
#include <unordered_map>

#include "bed_utils.hpp"

class HaplotypeShifter {

public:

    void load_tracks() ;
    std::string shift_coordinate(std::string, int) ;

    std::unordered_map<std::string, std::vector<Track>> tracks[2] ;
    std::unordered_map<std::string, std::vector<std::pair<int, int>>> offsets ;

    Track* find(std::string chrom, int position, int len) ;

} ;

#endif
