#include <omp.h>
#include <iomanip>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include <string_view>

#include "aggregator.hpp"
#include "chromosomes.hpp"

using namespace std ;

void Aggregator::run() {
    auto c = Configuration::getInstance() ;
    load_sequences() ;
    dump_sequences() ;
}

void Aggregator::load_sequences() {
    auto c = Configuration::getInstance() ;
    num_batches = c->aggregate_batches + 1 ;
    cout << "Loading sequences from " << num_batches << " batches.." << endl ;
    vector<unordered_map<string, int>> _sequences(num_batches) ;
    int e = 0 ;
    #pragma omp parallel for num_threads(num_batches)
    for (int j = 0; j <= num_batches; j++) {
        string s_j = std::to_string(j) ;
        string path = c->workdir + "/solution_batch_" + s_j + ".fastq" ;
        ifstream txt_file(path) ;
        int i = 0 ;
        int count ;
        string line ;
        string name ;
        while (std::getline(txt_file, line)) {
            if (i == 0) {
                int p = line.rfind(':') ;
                name = line.substr(0, p - 1) ;
                count = std::stoi(line.substr(p + 1, line.length() - (p + 1))) ;
            }
            if (i == 1) {
                //cout << name << " " << count << endl ;
                string canon = canonicalize(line) ;
                if (_sequences[j].find(canon) == _sequences[j].end()) {
                    _sequences[j][canon] == 0 ;
                }
                _sequences[j][canon] += count ;
            }
            i++ ;
            i %= 4 ;
        }
        txt_file.close() ;
    }
    cout << "Merging.." << endl ;
    for (int j = 0; j <= num_batches; j++) {
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
    int i = 0 ;
    ofstream o(path) ;
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
