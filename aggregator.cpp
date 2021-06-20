#include <omp.h>
#include <iomanip>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <assert.h>

#include "aggregator.hpp"

using namespace std ;

void Aggregator::run() {
    auto c = Configuration::getInstance() ;
    find_high_abundance_sequences() ;
    load_sequences() ;
    dump_sequences() ;
}

void Aggregator::find_high_abundance_sequences() {
    auto c = Configuration::getInstance() ;
    num_batches = c->aggregate_batches ;
    cout << "Finding high-abundance strings from " << num_batches << " batches.." << endl ;
    vector<unordered_map<int, int>> _sequences(num_batches) ;
    int e = 0 ;
    // cout first pass
    #pragma omp parallel for num_threads(num_batches)
    for (int j = 0; j < num_batches; j++) {
        string s_j = std::to_string(j) ;
        string path = c->workdir + "/solution_batch_" + s_j + ".fastq" ;
        ifstream fastq_file(path) ;
        int i = 0 ;
        int count ;
        string line ;
        string name ;
        while (std::getline(fastq_file, line)) {
            if (i == 0) {
                int p = line.rfind('#') ;
                count = std::stoi(line.substr(p + 1, line.length() - (p + 1))) ;
            }
            if (i == 1) {
                string canon = canonicalize(line) ;
                int hash = std::hash<std::string>()(canon) ;
                if (_sequences[j].find(hash) == _sequences[j].end()) {
                    _sequences[j][hash] = 0 ;
                }
                _sequences[j][hash] += count ;
            }
            i++ ;
            i %= 4 ;
        }
        fastq_file.close() ;
    }
    cout << "Merging batches.." << endl ;
    for (int j = 0; j < num_batches; j++) {
        cout << "Batch " << j << " with " << _sequences[j].size() << " sequences." << endl ;
        for (auto it = _sequences[j].begin(); it != _sequences[j].end(); it++) {
            if (sequence_index.find(it->first) == sequence_index.end()) {
                sequence_index[it->first] = 0 ;
            }
            sequence_index[it->first] += it->second ;
        }
        _sequences[j].clear() ;
    }
    cout << "Loaded " << sequence_index.size() << " SFS sequences." << endl ;
    auto it = sequence_index.begin();
    while (it != sequence_index.end()) {
        if (it->second < c->cutoff) {
            it = sequence_index.erase(it) ; 
        } else {
            it++ ;
        }
    }
    cout << sequence_index.size() << " high-abundance sequences found." << endl ;
}

void Aggregator::load_sequences() {
    auto c = Configuration::getInstance() ;
    num_batches = c->aggregate_batches ;
    cout << "Loading sequences from " << num_batches << " batches.." << endl ;
    vector<unordered_map<string, int>> _sequences(num_batches) ;
    int e = 0 ;
    // cout first pass
    #pragma omp parallel for num_threads(num_batches)
    for (int j = 0; j < num_batches; j++) {
        string s_j = std::to_string(j) ;
        string path = c->workdir + "/solution_batch_" + s_j + ".fastq" ;
        ifstream fastq_file(path) ;
        int i = 0 ;
        int count ;
        string line ;
        while (std::getline(fastq_file, line)) {
            if (i == 0) {
                int p = line.rfind('#') ;
                count = std::stoi(line.substr(p + 1, line.length() - (p + 1))) ;
            }
            if (i == 1) {
                string canon = canonicalize(line) ;
                int hash = std::hash<std::string>()(canon) ;
                if (sequence_index.find(hash) != sequence_index.end()) {
                    if (_sequences[j].find(canon) == _sequences[j].end()) {
                        _sequences[j][canon] = 0 ;
                    }
                    _sequences[j][canon] += count ;
                }
            }
            i++ ;
            i %= 4 ;
        }
        fastq_file.close() ;
    }
    cout << "Merging.." << endl ;
    for (int j = 0; j < num_batches; j++) {
        cout << "Batch " << j << " with " << _sequences[j].size() << " sequences." << endl ;
        for (auto it = _sequences[j].begin(); it != _sequences[j].end(); it++) {
            if (sequences.find(it->first) == sequences.end()) {
                sequences[it->first] = 0 ;
            }
            sequences[it->first] += it->second ;
        }
        _sequences[j].clear() ;
    }
    cout << "Loaded " << sequences.size() << " seuqences." << endl ;
}

void Aggregator::dump_sequences() {
    auto c = Configuration::getInstance() ;
    cout << "Dumping aggregated counts.." << endl ;
    string path = c->workdir + "/solution_aggregated.fastq" ;
    ofstream o(path) ;
    int i = 0 ;
    int n = 0 ;
    for (auto it = sequences.begin(); it != sequences.end(); it++) {
        if (it->second >= c->cutoff) {
            o << "@sol_" << i << "#" << it->second << endl ;
            o << it->first << endl ;
            o << "+" << endl ;
            for (int j = 0; j < it->first.length(); j++) {
                o << "I" ;
            }
            o << endl ;
            n += 1 ;
        }
        i += 1 ;
    }
    cout << n << " useful sequences remaining." << endl ;
}

