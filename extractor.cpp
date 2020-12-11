#include <omp.h>
#include <iomanip>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <bits/stdc++.h> 

#include "config.hpp"
#include "extractor.hpp"
#include "chromosomes.hpp"

using namespace std ;

void Extractor::run() {
    cout << "Extractor primed.." << endl ;
    load_strings() ;
    load_read_ids() ;
    load_reads() ;
    //dump_reads() ;
}

void Extractor::load_strings() {
    auto c = Configuration::getInstance() ;
    // load strings
    ifstream txt_file(c->bed) ;
    string line ;
    int n = 0 ;
    string name ;
    while (std::getline(txt_file, line)) {
        if (n == 0) {
            name = line.substr(1, line.size() - 1) ;
        }
        if (n == 1) {
            sequences[line].name = name ;
        }
        n += 1 ;
        n %= 4 ;
    }
    txt_file.close() ;
    cout << "Loaded " << sequences.size() << " strings.." << endl ;
}

void Extractor::load_read_ids() {
    auto c = Configuration::getInstance() ;
    int l = 0 ;
    int n = 0 ;
    string seq ;
    string line ;
    ifstream txt_file(c->fasta) ;
    while (std::getline(txt_file, line)) {
        if (n == 0) {
            //seq = line.substr(1, line.size() - 1) ;
            seq = line ;
        }
        if (n == 1) {
            if (sequences.find(seq) != sequences.end()) {
                stringstream ss(line) ;
                string token ;
                while(getline(ss, token, '$')) {
                    sequences[seq].read_ids.push_back(token) ;
                    read_ids[token].push_back(seq) ;
                    l += 1 ;
                }
            }
        }
        n += 1 ;
        n %= 2 ;
    }
    txt_file.close() ;
    cout << "Loaded " << read_ids.size() << " read ids.." << endl ;
}

void Extractor::load_reads() {
    auto c = Configuration::getInstance() ;
    gzFile fastq_file ;
    kseq_t* fastq_iterator ;
    fastq_file = gzopen(c->fastq.c_str(), "r") ;
    fastq_iterator = kseq_init(fastq_file) ;
    int l = 0 ;
    int n = 0 ;
    ofstream out_fastq_file(c->workdir + "/selected_reads.fastq") ;
    while ((l = kseq_read(fastq_iterator)) >= 0) {
        auto f = fastq_entry_t(fastq_iterator->name.s, fastq_iterator->seq.s, fastq_iterator->qual.s) ;
        if (read_ids.find(f.head) != read_ids.end()) {
            //reads[f.head] = f ;
            //for (auto s: read_ids[fastq_iterator->name.s]) {
            //    auto ss = strstr(f.seq.c_str(), s.c_str()) ;
            //    auto sr = strstr(f.seq.c_str(), reverse_complement(s).c_str()) ;
            //    if (ss != nullptr || sr != nullptr) { 
            //        reads[s].push_back(f) ;
            //        n += 1 ;
            //    }
            //}
            out_fastq_file << "@" + f.head << endl ;
            out_fastq_file << f.seq << endl ;
            out_fastq_file << "+" << endl ;
            out_fastq_file << f.qual << endl ;
            n += 1 ;
        }
        //if (n == 100) {
        //    break ;
        //}
    }
    cout << "Loaded " << n << " reads.. " << endl ;
}

//void Extractor::dump_reads() {
//    auto c = Configuration::getInstance() ;
//    //#pragma omp parallel for num_threads(48)
//    for (size_t bucket = 0; bucket < reads.bucket_count(); bucket++) {
//        for (auto it = reads.begin(bucket); it != reads.end(bucket); it++) {
//            ofstream fastq_file(c->workdir + "/seqs/" + std::to_string(int(bucket)) + ".fastq") ;
//            for (const auto fastq_entry : it->second) {
//                fastq_file << "@" << fastq_entry.head << endl 
//                    << fastq_entry.seq << endl
//                    << "+" << endl
//                    << fastq_entry.qual << endl ;
//            }
//        }
//    }
//}
