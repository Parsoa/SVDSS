#include <omp.h>
#include <iomanip>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <bits/stdc++.h> 

#include "config.hpp"
#include "chromosomes.hpp"
#include "extractor.hpp"

using namespace std ;

void Extractor::run() {
    cout << "Extractor primed.." << endl ;
    load_strings() ;
    load_read_ids() ;
    load_reads() ;
    dump_reads() ;
}

void Extractor::load_strings() {
    auto c = Configuration::getInstance() ;
    // load strings
    ifstream txt_file(c->workdir + "/short_unmapped.bed") ;
    string line ;
    while (std::getline(txt_file, line)) {
        stringstream ss(line) ;
        vector<string> tokens ;
        string token ;
        while (getline(ss, token, '\t')) {
            tokens.push_back(token) ;
        }
        sequences[tokens[3]] = Locus{tokens[1], std::stoi(tokens[2]), std::stoi(tokens[4])} ;
    }
    txt_file.close() ;
    cout << "Loaded " << sequences.size() << " strings.." << endl ;
}

void Extractor::load_read_ids() {
    auto c = Configuration::getInstance() ;
    gzFile fastq_file ;
    kseq_t* fastq_iterator ;
    string path = c->workdir + "/solution.fastq" ;
    fastq_file = gzopen(path.c_str(), "r") ;
    fastq_iterator = kseq_init(fastq_file) ;
    int l = 0 ;
    int n = 0 ;
    while ((l = kseq_read(fastq_iterator)) >= 0) {
        if (sequences.find(fastq_iterator->seq.s) != sequences.end() || sequences.find(reverse_complement(fastq_iterator->seq.s)) != sequences.end()) {
            stringstream ss(fastq_iterator->name.s) ;
            vector<string> tokens ;
            string token ;
            while (getline(ss, token, '.')) {
                tokens.push_back(token) ;
            }
            read_ids[tokens[0]].push_back(fastq_iterator->seq.s) ;
            //read_ids[tokens[0]].push_back(reverse_complement(fastq_iterator->seq.s)) ;
            n += 1 ;
        }
    }
    cout << "Loaded " << n << " read ids.." << endl ;
}

void Extractor::load_reads() {
    auto c = Configuration::getInstance() ;
    gzFile fastq_file ;
    kseq_t* fastq_iterator ;
    string path = c->workdir + "/reads.fastq" ;
    fastq_file = gzopen(path.c_str(), "r") ;
    fastq_iterator = kseq_init(fastq_file) ;
    int l = 0 ;
    int n = 0 ;
    while ((l = kseq_read(fastq_iterator)) >= 0) {
        auto f = fastq_entry_t(fastq_iterator->name.s, fastq_iterator->seq.s, fastq_iterator->qual.s) ;
        if (read_ids.find(fastq_iterator->name.s) != read_ids.end()) {
            for (auto s: read_ids[fastq_iterator->name.s]) {
                reads[s].push_back(f) ;
                n += 1 ;
            }
            //cout << n << endl ;
        }
        //if (n == 100) {
        //    break ;
        //}
    }
    cout << "Loaded " << n << " reads.. " << endl ;
}

void Extractor::dump_reads() {
    auto c = Configuration::getInstance() ;
    //#pragma omp parallel for num_threads(48)
    for (size_t bucket = 0; bucket < reads.bucket_count(); bucket++) {
        for (auto it = reads.begin(bucket); it != reads.end(bucket); it++) {
            ofstream fastq_file(c->workdir + "/seqs/" + std::to_string(int(bucket)) + ".fastq") ;
            for (const auto fastq_entry : it->second) {
                fastq_file << "@" << fastq_entry.head << endl 
                    << fastq_entry.seq << endl
                    << "+" << endl
                    << fastq_entry.qual << endl ;
            }
        }
    }
}
