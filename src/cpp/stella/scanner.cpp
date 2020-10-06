#include <omp.h>
#include <iomanip>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include <string_view>

#include "shifter.hpp"
#include "scanner.hpp"
#include "chromosomes.hpp"

using namespace std ;

void Scanner::run() {
    auto c = Configuration::getInstance() ;
    load_chromosomes(c->reference) ;
    load_sequences() ;
    scan_reference(c->threads) ;
    dump_sequences() ;
}

void Scanner::load_sequences() {
    auto c = Configuration::getInstance() ;
    cout << "Loading sequences.." << endl ;
    vector<unordered_map<string, int>> _sequences(c->batch_size) ;
    #pragma omp parallel for num_threads(48)
    for (int j = 0; j < c->batch_size; j++) {
        string s_j = std::to_string(j) ;
        string index = j < 10 ? "00" + s_j : j < 100 ? "0" + s_j : s_j ;
        string path = c->workdir + "/solution/solution.mmp.unmapped.batch_" + index + ".txt" ;
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
            string name = tokens[0].substr(0, p - 1) ;
            int count = std::stoi(tokens[0].substr(p + 1, tokens[0].length() - (p + 1))) ;
            if (count > 5) {
                if (tokens[3].length() >= 10 && tokens[3].length() <= 40) {
                    _sequences[j][tokens[3]] = count ;
                }
            }
            //cout << "#" << tokens[3] << "#" ;
        }
        txt_file.close() ;
    }
    cout << "Merging.." << endl ;
    for (int j = 0; j < c->batch_size; j++) {
        //cout << "Batch " << j << " with " << _sequences[j].size() << " sequences." << endl ;
        sequences.insert(_sequences[j].begin(), _sequences[j].end()) ;
        _sequences[j].clear() ;
    }
    cout << "Loaded " << sequences.size() << " seuqences." << endl ;
}

void Scanner::scan_reference(int threads) {
    cout << "--------------------------------------------------------- " << endl ;
    cout << "Scanning reference genome.." << endl ;
    threads = min(threads, int(chromosome_seqs.size())) ;
    for (int t = 0; t < threads - 1; t++) {
        cout << endl ;
    }
    int m = 0 ;
    int n = chromosomes.size() ;
    vector<batch_type> batches(chromosomes.size()) ; 
    while (m < n) {
        int p = m ;
        #pragma omp parallel for num_threads(threads)
        for (int i = 0; i < threads; i++) {
            if (m + i < n) {
                batches[i + m] = scan_chromosome(chromosomes[m + i], threads, omp_get_thread_num()) ;
            }
        }
        m += threads ;
        if (m >= chromosomes.size()) {
            m = n ;
        }
    }
    cout << "---------------------------------------------------------" << endl ;
    for (int i = 0; i < chromosomes.size(); i++) {
        cout << "Merging " << std::left << std::setw(6) << chromosomes[i] << ", " << std::setw(9) << batches[i].size() << " matches.." << endl ;
        for(auto it = batches[i].begin(); it != batches[i].end(); it++) {
            kmers[it->first].insert(kmers[it->first].end(), it->second.begin(), it->second.end()) ;
        }
        batches[i].clear() ;
    }
    cout << "Reference scan completed." << endl ;
}

batch_type Scanner::scan_chromosome(string chrom, int threads, int index) {
    batch_type _kmers ;
    char* seq = chromosome_seqs[chrom] ;
    char* kmer = (char*) malloc(sizeof(char) * 41) ;
    uint64_t l = strlen(chromosome_seqs[chrom]) ;
    for (int i = 2000; i < l - 2000; i++) {
        for (int k = 10; k < 40; k++) {
            strncpy(kmer, seq + i, k) ;
            kmer[k + 1] = '\0' ;
            if (sequences.find(kmer) != sequences.end()) {
                _kmers[kmer].push_back(Locus{chrom, i}) ;
            }
        }
        if (i % 10000000 == 0) {
            cout_mutex.lock() ;
            for (int j = 0; j < (threads - 1) - index; j++) {
                cout << "\x1b[A" ;
            }
            cout << "\r"<< std::left << std::setw(6) << chrom << " progress " << std::fixed << std::setprecision(3) << float(i) / float(l) ;
            for (int j = 0; j < (threads - 1) - index; j++) {
                cout << endl ;
            }
            cout_mutex.unlock() ;
        }
    }
    return _kmers ;
}

void Scanner::dump_sequences() {
    auto c = Configuration::getInstance() ;
    auto shifter = Shifter() ;
    shifter.load_tracks() ;
    ofstream seq_file(c->workdir + "/unmapped_seqs.bed") ;
    ofstream track_file(c->workdir + "/unmapped_tracks.bed") ;
    ofstream coord_file(c->workdir + "/unmapped_coordinates.bed") ;
    unordered_map<Track, int> unmapped_tracks ;
    for (auto it = kmers.begin(); it != kmers.end(); it++) {
        for (auto& locus: it->second) {
            auto match = shifter.find(locus.chrom, locus.position, it->first.length()) ;
            if (match != nullptr) {
                Track& track = *match ;
                int count = sequences[it->first] ;
                if (unmapped_tracks.find(track) == unmapped_tracks.end()) {
                    unmapped_tracks[track] = 0 ;
                    track_file << track.chrom << "\t" << track.begin << "\t" << track.end << "\t" << it->first << "\t" << count << endl ;
                }
                seq_file << it->first << "\t" << count << "\t" << locus.chrom << "\t" << locus.position << "\t" << track.get_name() << endl ;
            } else {
                auto pos = shifter.shift_coordinate(locus.chrom, locus.position) ;
                coord_file << locus.chrom << "\t" << locus.position << "\t" << pos << "\t" << locus.count << "\t" << it->first << endl ;
            }
        }
    }
}
