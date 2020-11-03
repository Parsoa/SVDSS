#include <omp.h>
#include <iomanip>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <bits/stdc++.h> 

#include "config.hpp"
#include "shifter.hpp"
#include "intersector.hpp"

using namespace std ;

void Intersector::run() {
    auto c = Configuration::getInstance() ;
    cout << c->workdir << endl ;
    auto sequences = load_sequences(c->workdir) ;
    auto corrected_sequences = load_sequences(c->workdir + "-edited") ;
    for (auto sequence: sequences) {
        if (corrected_sequences.find(sequence.first) != corrected_sequences.end()) {
            shared_sequences[sequence.first] = sequence.second ; 
        }
    }
    cout << shared_sequences.size() << " common sequences." << endl ;
    dump_sequences() ;
}

unordered_map<string, Locus> Intersector::load_sequences(string base) {
    auto c = Configuration::getInstance() ;
    cout << "Loading sequences.." << endl ;
    unordered_map<string, Locus> sequences ;
    vector<unordered_map<string, Locus>> _sequences(c->batch_size) ;
    #pragma omp parallel for num_threads(48)
    for (int j = 0; j < c->batch_size; j++) {
        string s_j = std::to_string(j) ;
        string index = j < 10 ? "00" + s_j : j < 100 ? "0" + s_j : s_j ;
        string path = base + "/" + c->type + "/solution.mapped.batch_" + index + ".txt" ;
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
    return sequences ;
}

void Intersector::dump_sequences() {
    auto c = Configuration::getInstance() ;
    auto shifter = Shifter() ;
    shifter.load_tracks() ;

    ofstream mapped_bed(c->workdir + "/" + c->type + "_mapped.intersect.bed") ;
    ofstream mapped_fasta(c->workdir + "/" + c->type + "_mapped.intersect.fasta") ;
    ofstream mapped_bed_seqs(c->workdir + "/" + c->type + "_mapped.seqs.intersect.bed") ;

    ofstream unmapped_bed(c->workdir + "/" + c->type + "_unmapped.intersect.bed") ;
    ofstream unmapped_fasta(c->workdir + "/" + c->type + "_unmapped.intersect.fasta") ;
    ofstream unmapped_fastq(c->workdir + "/" + c->type + "_unmapped.intersect.fastq") ;

    unordered_map<Track, int> mapped_tracks ;
    int m = 0 ;
    for (auto it = shared_sequences.begin(); it != shared_sequences.end(); it++) {
        auto& locus = it->second ;
        auto match = shifter.find(locus.chrom, locus.position, it->first.length()) ;
        if (match != nullptr) {
            auto& track = *match ;
            auto pos = shifter.shift_coordinate(locus.chrom, locus.position) ;
            if (mapped_tracks.find(track) == mapped_tracks.end()) {
                m += 1 ;
                mapped_tracks[track] = 0 ;
                mapped_bed << track.chrom << "\t" << track.begin << "\t" << track.end << "\t" << it->first << "\t" << locus.count << "\t" << locus.position << "\t" << pos << endl ;
            }
            mapped_fasta << ">" << locus.chrom << "_" << locus.position << endl ; 
            mapped_fasta << it->first << endl ;
            mapped_bed_seqs << track.chrom << "\t" << track.begin << "\t" << track.end << "\t" << it->first << "\t" << locus.count << "\t" << locus.position << endl ;
        } else {
            auto pos = shifter.shift_coordinate(locus.chrom, locus.position) ;
            unmapped_bed << locus.chrom << "\t" << locus.position << "\t" << pos << "\t" << it->first << "\t" << locus.count << endl ;
            // FASTQ file
            unmapped_fastq << "@" << locus.chrom << "_" << locus.position << endl ;
            unmapped_fastq << it->first << endl ;
            unmapped_fastq << "+" << locus.chrom << "_" << locus.position << endl ;
            for (int j = 0; j < it->first.length(); j++) {
                unmapped_fastq << "I" ;
            }
            unmapped_fastq << endl ;
            // FASTA file
            unmapped_fasta << ">" << locus.chrom << "_" << locus.position << endl ;
            unmapped_fasta << it->first << endl ;
        }
    }
}
