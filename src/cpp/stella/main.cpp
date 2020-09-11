#include <string>
#include <cstdlib>
#include <cstdint>
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
#include "aggregator.hpp"
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
    auto c = Configuration::getInstance() ;
    if (strcmp(argv[1], "search") == 0) {
        c->parse(argc - 2, argv + 2) ;
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
    c->parse(argc - 1, argv + 1) ;
    create_workdir() ;
    if (strcmp(argv[1], "scan") == 0) {
        auto scanner = new Scanner() ;
        scanner->run() ;
    }
    if (strcmp(argv[1], "find") == 0) {
        auto finder = new Finder() ;
        finder->run() ;
    }
    if (strcmp(argv[1], "aggregate") == 0) {
        auto aggregator = new Aggregator() ;
        aggregator->run() ;
    }
}

