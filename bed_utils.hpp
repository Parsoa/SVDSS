#ifndef BED_HPP
#define BED_HPP

#include <vector>
#include <string>
#include <unordered_map>
#include <sstream>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <fstream>
#include <stdlib.h>
#include <assert.h>

#include "lprint.hpp"

using namespace std ;

struct Locus {
    std::string chrom ;
    int position ;
    int count ;
} ;

struct Track {
    std::string chrom ;
    uint64_t begin ;
    uint64_t end ;
    std::string svtype ;
    int svlen ;

    bool operator==(const Track &o) const {
        return chrom == o.chrom && begin == o.begin && end == o.end && svtype == o.svtype ;
    }

    bool operator!=(const Track &o) const {
        return !(*this == o) ;
    }

    std::string get_name() const ;
};

namespace std {
    template <> struct hash<Track> {
        std::size_t operator()(const Track& t) const {
            uint64_t hash = 0 ;
            hash += t.begin ;
            hash += t.end << 32 ;
            return std::hash<uint64_t>()(hash) ;
        }
    };
}

bool track_comparator(Track a, Track b) ;

std::vector<Track> load_tracks_from_file(std::string path) ;
std::unordered_map<Track, int> load_tracks_from_file_as_dict(std::string path) ;

#endif
