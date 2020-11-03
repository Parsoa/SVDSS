#include <string>
#include <fstream>
#include <stdlib.h>
#include <iostream>
#include <assert.h>
#include <unordered_map>

#include "bed_utils.hpp"

using namespace std ;

bool track_comparator(Track a, Track b) {
    if (a.chrom[3] == b.chrom[3]) {
        return a.begin < b.end ;
    } else {
        return a.chrom[3] < b.chrom[3] ;
    }
}

string Track::get_name() const {
    return svtype + "@" + chrom + "_" + std::to_string(begin) + "_" + std::to_string(end) ;
}

std::vector<Track> load_tracks_from_file(string path) {
    cout << "Parsing BED file: " << path << ".." << endl ;
    std::ifstream bed_file(path) ;
    std::string line ;
    int i = 0 ;
    std::vector<Track> tracks ;
    unordered_map<string, int> header ;
    while (std::getline(bed_file, line)) {
        istringstream iss(line) ;
        if (i == 0) {
            if (line[0] != '#') {
                cout << "BED header not present (first line doesn't begin with #). Aborting.." << endl ;
            }
            vector<string> tokens{istream_iterator<string>{iss}, istream_iterator<string>{}} ;
            for (int i = 0; i < tokens.size(); i++) {
                header[tokens[i]] = i ;
            }
            i += 1 ;
        } else {
            vector<string> tokens{istream_iterator<string>{iss}, istream_iterator<string>{}} ;
            Track track ;
            track.chrom = tokens[0] ;
            track.begin = std::stoi(tokens[1]) ;
            track.end = std::stoi(tokens[2]) ;
            track.svtype = tokens[3] ;
            track.svlen = std::stoi(tokens[4]) ;
            tracks.push_back(track) ;
        }
    }
    cout << "Loaded " << tracks.size() << " tracks." << endl ;
    return tracks ;
}

std::unordered_map<Track, int> load_tracks_from_file_as_dict(string path) {
    auto _tracks = load_tracks_from_file(path) ;
    std::unordered_map<Track, int> tracks ;
    for (auto track = _tracks.begin(); track != _tracks.end(); track++) {
        tracks[*track] = 1 ;
    }
    return tracks ;
}

