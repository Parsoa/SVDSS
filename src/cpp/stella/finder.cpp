#include <omp.h>
#include <iomanip>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <bits/stdc++.h> 

#include "config.hpp"
#include "finder.hpp"
#include "shifter.hpp"

using namespace std ;

void Finder::run() {
    auto c = Configuration::getInstance() ;
    load_sequences() ;
    dump_sequences() ;
}

void Finder::load_sequences() {
    auto c = Configuration::getInstance() ;
    cout << "Loading sequences.." << endl ;
    vector<unordered_map<string, Locus>> _sequences(c->batch_size) ;
    #pragma omp parallel for num_threads(48)
    for (int j = 0; j < c->batch_size; j++) {
        string s_j = std::to_string(j) ;
        string index = j < 10 ? "00" + s_j : j < 100 ? "0" + s_j : s_j ;
        string path = c->workdir + "/solution/solution.mmp.mapped.batch_" + index + ".txt" ;
        ifstream txt_file(path) ;
        string line ;
        while (std::getline(txt_file, line)) {
            stringstream ss(line) ;
            vector<string> tokens ;
            string token ;
            while (getline(ss, token, ' ')) {
                tokens.push_back(token) ;
            }
            assert(tokens.size() == 4) ;
            int p = tokens[0].rfind('#') ;
            int count = std::stoi(tokens[0].substr(p + 1, tokens[0].length() - (p + 1))) ;
            if (count >= 5) {
                _sequences[j][tokens[3]] = Locus{tokens[1], std::stoi(tokens[2]), count} ;
            }
        }
        txt_file.close() ;
    }
    cout << "Merging.." << endl ;
    for (int j = 0; j < c->batch_size; j++) {
        sequences.insert(_sequences[j].begin(), _sequences[j].end()) ;
        _sequences[j].clear() ;
    }
    cout << "Loaded " << sequences.size() << " seuqences." << endl ;
}

void Finder::dump_sequences() {
    auto c = Configuration::getInstance() ;
    auto shifter = Shifter() ;
    shifter.load_tracks() ;
    ofstream seq_file(c->workdir + "/mapped_seqs.bed") ;
    ofstream track_file(c->workdir + "/mapped_tracks.bed") ;
    ofstream coord_file(c->workdir + "/mapped_coordinates.bed") ;
    ofstream fastq_file(c->workdir + "/mapped_coordinates.fastq") ;
    unordered_map<Track, int> mapped_tracks ;
    for (auto it = sequences.begin(); it != sequences.end(); it++) {
        auto& locus = it->second ;
        auto match = shifter.find(locus.chrom, locus.position, it->first.length()) ;
        if (match != nullptr) {
            auto& track = *match ;
            if (mapped_tracks.find(track) == mapped_tracks.end()) {
                mapped_tracks[track] = 0 ;
                track_file << track.chrom << "\t" << track.begin << "\t" << track.end << "\t" << track.svtype << "\t" << track.svlen << "\t" << it->first << "\t" << locus.count << endl ;
            }
            seq_file << it->first << "\t" << locus.count << "\t" << locus.chrom << "\t" << locus.position << "\t" << track.get_name() << endl ;
        } else {
            auto pos = shifter.shift_coordinate(locus.chrom, locus.position) ;
            coord_file << locus.chrom << "\t" << locus.position << "\t" << pos << "\t" << locus.count << "\t" << it->first << endl ;
            fastq_file << "@" << locus.chrom << "_" << locus.position << endl ;
            fastq_file << it->first << endl ;
            fastq_file << "+" << locus.chrom << "_" << locus.position << endl ;
            for (int j = 0; j < it->first.length(); j++) {
                fastq_file << "I" ;
            }
            fastq_file << endl ;
        }
    }
}
