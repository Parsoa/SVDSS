#include <string>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <iostream>
#include <algorithm>

#include "config.hpp"
#include "intersector.hpp"

using namespace std ;

void Intersector::load_tracks() {
    auto c = Configuration::getInstance() ;
    cout << "Intersector loading tracks.." << endl ;
    string path = c->bed ;
    std::ifstream bed_file(path) ;
    std::string line ;
    unordered_map<string, int> header ;
    while (std::getline(bed_file, line)) {
        if (line[0] == '#') {
            continue ;
        }
        istringstream iss(line) ;
        vector<string> tokens{istream_iterator<string>{iss}, istream_iterator<string>{}} ;
        Track track ;
        track.chrom = tokens[0] ;
        track.begin = std::stoi(tokens[1]) ;
        track.end = std::stoi(tokens[2]) ;
        track.svtype = tokens[3] ;
        if (track.svtype == "DEL") {
            track.svlen = track.end - track.begin ;
        } else if (track.svtype == "INS") {
            track.svlen = tokens[5].length() ;
        } else {
            track.svlen = abs(std::stoi(tokens[4])) ;
        }
        tracks[tokens[0]].push_back(track) ;
    }
    cout << "Loaded " << tracks.size() << " tracks." << endl ;
}

Track* Intersector::find(string chrom, int position, int len) {
    for (auto track = tracks[chrom].begin(); track != tracks[chrom].end(); track++) {
        if (chrom == track->chrom) {
            if (track->svtype == "INS") {
                if (position <= track->begin + track->svlen && position + len >= track->begin) {
                    return &(*track) ;
                }
            }
            if (track->svtype == "DEL") {
                if (position <= track->begin && position + len > track->begin) {
                    return &(*track) ;
                }
            }
            if (track->svtype == "INV") {
                if (position <= track->begin + track->svlen && position + len >= track->begin) {
                    return &(*track) ;
                }
            }
        }
    }
    return nullptr ;
}
