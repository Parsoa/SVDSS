#include <omp.h>
#include <string>
#include <cstdlib>
#include <cstdint>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <assert.h>
#include <iostream>
#include <stdlib.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unordered_map>

#include "config.hpp"
#include "finder.hpp"
#include "shifter.hpp"
#include "scanner.hpp"
#include "extractor.hpp"
#include "ping_pong.hpp"
#include "aggregator.hpp"
#include "intersector.hpp"
#include "chromosomes.hpp"

using namespace std ;

// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\

void create_workdir() {
    auto c = Configuration::getInstance() ;
    struct stat info ;
    if (stat(c->workdir.c_str(), &info) != 0) {
        cout << "Working directory " << c->workdir << " does not exist. Aborting.." << endl ;
        exit(0) ;
    }
}

int main(int argc, char** argv) {
    cout << "Stella, mapping-free variation discovery." << endl ;
    if (strcmp(argv[1], "pingpong") == 0) {
        auto pingpong = new PingPong() ;
        int retcode = pingpong->run(argc - 1, argv + 1) ;
        exit(retcode) ;
    }
    auto c = Configuration::getInstance() ;
    c->parse(argc - 1, argv + 1) ;
    if (strcmp(argv[1], "search") == 0) {
        load_chromosomes(c->reference) ;
        char* s = nullptr ;
        for (auto chrom = chromosome_seqs.begin(); chrom != chromosome_seqs.end(); chrom++) {
            s = strstr(chromosome_seqs[chrom->first], argv[2]) ;
            if (s != nullptr) {
                cout << chrom->first << ": Found at " << s - chromosome_seqs[chrom->first] << endl ;
            }
        }
        exit(0) ;
    }
    if (strcmp(argv[1], "bedsearch") == 0) {
        load_chromosomes(c->reference) ;
        std::ifstream bed_file(c->bed) ;
        std::string line ;
        vector<string> lines ;
        while (std::getline(bed_file, line)) {
            istringstream iss(line) ;
            vector<string> tokens{istream_iterator<string>{iss}, istream_iterator<string>{}} ;
            lines.push_back(tokens[3]) ;
        }
        vector<int> n ;
        vector<int> m ;
        for(int i = 0; i < 48; i++) {
            n.push_back(0) ;
            m.push_back(0) ;
        }
        cout << "Searching for " << lines.size() << " sequences.." << endl ;
        #pragma omp parallel for num_threads(48)
        for (int i = 0; i < lines.size(); i++) {
            int t = omp_get_thread_num() ;
            char* s = nullptr ;
            for (auto chrom = chromosome_seqs.begin(); chrom != chromosome_seqs.end(); chrom++) {
                s = strstr(chromosome_seqs[chrom->first], lines[i].c_str()) ;
                if (s != nullptr) {
                    m[t] += 1 ;
                    cout << "Match " << lines[i] << endl ;
                    break ;
                }
            }
            n[t] += 1 ;
        }
        int _n ;
        int _m ;
        for(int i = 0; i < 48; i++) {
            _n += n[i] ;
            _m += m[i] ;
        }
        cout << "Found " << _m << " out of " << _n << endl ;
        exit(0) ;
    }
    create_workdir() ;
    if (strcmp(argv[1], "scan") == 0) {
        auto scanner = new Scanner() ;
        scanner->run() ;
    }
    if (strcmp(argv[1], "find") == 0) {
        auto finder = new Finder() ;
        finder->run() ;
    }
    if (strcmp(argv[1], "extract") == 0) {
        auto extractor = new Extractor() ;
        extractor->run() ;
    }
    if (strcmp(argv[1], "aggregate") == 0) {
        auto aggregator = new Aggregator() ;
        aggregator->run() ;
    }
    if (strcmp(argv[1], "intersect") == 0) {
        auto intersector = new Intersector() ;
        intersector->run() ;
    }
}

