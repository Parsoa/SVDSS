#include <string>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <iostream>
#include <algorithm>

#include "config.hpp"
#include "shifter.hpp"
#include "haplotype_shifter.hpp"

using namespace std ;

void HaplotypeShifter::load_tracks() {
    auto c = Configuration::getInstance() ;
    cout << "HaplotypeShifter loading tracks.." << endl ;
    string path = c->vcf ;
    std::ifstream vcf_file(path) ;
    std::string line ;
    unordered_map<string, int> header ;
    vector<ofstream> out_files ;
    for (int i = 0; i < 2; i++) {
        out_files.emplace_back(ofstream{c->workdir + "/shifted_" + std::to_string(i + 1) + ".bed"}) ;
    }
    int pos ;
    int offset[] = {0, 0} ;
    string chrom = "null" ;
    while (std::getline(vcf_file, line)) {
        if (line[0] == '#') {
            continue ;
        }
        istringstream iss(line) ;
        vector<string> tokens{istream_iterator<string>{iss}, istream_iterator<string>{}} ;
        if (tokens[0] != chrom) {
            offset[0] = 0 ;
            offset[1] = 0 ;
            chrom = tokens[0] ;
            cout << chrom << endl ;
        }
        pos = std::stoi(tokens[1]) ;
        auto ref = tokens[3] ;
        istringstream alt_iss(tokens[4]) ;
        vector<string> alt_tokens ;
        string alt_token ;
        while (std::getline(alt_iss, alt_token, ',')) {
            alt_tokens.push_back(alt_token) ;
        }
        auto genotype = tokens[11] ; // HG00733
        //cout << tokens[4] << " " << alt_tokens.size() << " " << genotype << endl ;
        for (int i = 0; i < 2; i++) {
            Track track ;
            track.chrom = tokens[0] ;
            track.begin = pos + offset[i] ;
            if (genotype[i * 2] != '0') {
                int alt = std::stoi(genotype.substr(i * 2, 1)) - 1 ;
                offset[i] += alt_tokens[alt].length() - ref.length() ;
                track.end = track.begin + alt_tokens[alt].length() ; 
                track.svlen = alt_tokens[alt].length() - ref.length() ;
                tracks[i][chrom].push_back(track) ;
                out_files[i] << track.chrom << "\t" << track.begin << "\t" << track.end << "\t" << ref << "\t" << alt_tokens[alt] << "\t" << track.svlen << "\t" << genotype << "\t" << tokens[0] << "_" << tokens[1] << endl ;
            }
        }
    }
}

Track* HaplotypeShifter::find(string chrom, int position, int len) {
    //for (auto track = tracks[chrom].begin(); track != tracks[chrom].end(); track++) {
    //    if (position <= track->begin + track->svlen && position + len >= track->begin) {
    //        return &track->second ;
    //    }
    //}
    return nullptr ;
}

string HaplotypeShifter::shift_coordinate(string chrom, int pos) {
    int offset = 0 ;
    for (auto it: offsets[chrom]) {
        if (it.first > pos) {
            return std::to_string(pos - offset) ;
        }
        offset = it.second ;
    }
    return std::to_string(pos - offset) ;
}

