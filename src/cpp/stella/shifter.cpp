#include <string>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <iostream>
#include <algorithm>

#include "config.hpp"
#include "shifter.hpp"

using namespace std ;

void Shifter::load_tracks() {
    auto c = Configuration::getInstance() ;
    cout << "Shifter loading tracks.." << endl ;
    string path = c->bed ;
    std::ifstream bed_file(path) ;
    std::string line ;
    unordered_map<string, int> header ;
    vector<ofstream> out_files ;
    for (int i = 0; i < 2; i++) {
        out_files.emplace_back(ofstream{c->workdir + "/shifted_" + std::to_string(i + 1) + ".bed"}) ;
    }
    int offset[] = {0, 0} ;
    string chrom = "null" ;
    int end = 0 ;
    int begin = 0 ;
    string previous ;
    int n = 0 ;
    vector<string> chromosomes {"chr1", "chr2", "chr3", "chr4", "chr5"} ;
    for (auto chrom: chromosomes) {
        for (int i = 0; i < 2; i++) {
            offsets[chrom + "_" + std::to_string(i)].push_back(std::make_pair(0, 0)) ;
        }
    }
    while (std::getline(bed_file, line)) {
        if (line[0] == '#') {
            continue ;
        }
        istringstream iss(line) ;
        vector<string> tokens{istream_iterator<string>{iss}, istream_iterator<string>{}} ;
        Track track ;
        if (std::find(chromosomes.begin(), chromosomes.end(), tokens[0]) == chromosomes.end()) {
            continue ;
        }
        n += 1 ;
        if (tokens[0] != chrom) {
            offset[0] = 0 ;
            offset[1] = 0 ;
            chrom = tokens[0] ;
            //cout << chrom << endl ;
            end = 0 ;
        }
        begin = std::stoi(tokens[1]) ;
        if (begin < end) {
            //cout << "Overlap:" << endl ;
            //cout << line << endl ;
            //cout << previous << endl ;
        }
        end = std::stoi(tokens[2]) ;
        track.chrom = tokens[0] ;
        track.begin = begin ;
        track.end = end ;
        track.svtype = tokens[3] ;
        if (track.svtype == "DEL") {
            track.svlen = track.end - track.begin ;
        } else if (track.svtype == "INS") {
            track.svlen = tokens[5].length() ;
        } else {
            track.svlen = abs(std::stoi(tokens[4])) ;
        }
        for (int i = 0; i < 2; i++) {
            Track alt ;
            alt.chrom = tokens[0] + "_" + std::to_string(i + 1) ;
            alt.begin = track.begin + offset[i] ;
            alt.end = track.end + offset[i] ;
            alt.svtype = track.svtype ;
            alt.svlen = track.svlen ;
            if (tokens[10 + i] == "1") {
                if (tokens[3] == "DEL") {
                    offset[i] -= (track.end - track.begin) ;
                } else if (tokens[3] == "INS") {
                    offset[i] += (tokens[5].length() - 1) ; 
                } else {
                    // inversions don't change offset
                }
                tracks[alt] = track ;
                out_files[i] << alt.chrom << "\t" << alt.begin << "\t" << alt.end << "\t" << alt.svtype << "\t" << alt.svlen << endl ;
                offsets[alt.chrom].push_back(std::make_pair(alt.begin, offset[i])) ;
            }
        }
        previous = line ;
    }
    cout << "Loaded " << tracks.size() << " tracks." << endl ;
}

Track* Shifter::find(string chrom, int position, int len) {
    for (auto track = tracks.begin(); track != tracks.end(); track++) {
        //if (chrom.find(track->first.chrom) != string::npos) {
        if (chrom == track->first.chrom) {
            if (track->first.svtype == "INS") {
                if (position <= track->first.begin + track->first.svlen && position + len >= track->first.begin) {
                    return &track->second ;
                }
            }
            if (track->first.svtype == "DEL") {
                if (position <= track->first.begin && position + len > track->first.begin) {
                    return &track->second ;
                }
            }
        }
    }
    return nullptr ;
}

string Shifter::shift_coordinate(string chrom, int pos) {
    int i = 0 ;
    int offset = 0 ;
    for (auto it: offsets[chrom]) {
        if (it.first > pos) {
            //return pos + offset ;
            return std::to_string(pos - offset) ;
        }
        offset = it.second ;
    }
    //return pos + offset ;
    return std::to_string(pos - offset) ;
}

