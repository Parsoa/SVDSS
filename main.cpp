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
#include "ping_pong.hpp"
#include "extractor.hpp"
#include "aggregator.hpp"
#include "chromosomes.hpp"

using namespace std ;
//using namespace fs = std::filesystem ;

// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\

void create_workdir() {
    auto c = Configuration::getInstance() ;
    struct stat info ;
    cout << c->workdir << endl ;
    if (stat(c->workdir.c_str(), &info) != 0) {
        cout << "Working directory does not exist. creating.." << endl ;
        int check = mkdir(c->workdir.c_str(), 0777) ;
        if (check != 0) {
            cerr << "Failed to create output directory, aborting.." << endl ;
            exit(check) ;
        }
    }
}

int main(int argc, char** argv) {
    cout << "Ping-pong, comparative genome analysis using sample-specific string detection in accurate long reads." << endl ;
    auto c = Configuration::getInstance() ;
    //if (strcmp(argv[1], "search") == 0) {
    //    load_chromosomes(c->reference) ;
    //    char* s = nullptr ;
    //    for (auto chrom = chromosome_seqs.begin(); chrom != chromosome_seqs.end(); chrom++) {
    //        s = strstr(chromosome_seqs[chrom->first], argv[2]) ;
    //        if (s != nullptr) {
    //            cout << chrom->first << ": Found at " << s - chromosome_seqs[chrom->first] << endl ;
    //        }
    //    }
    //    exit(0) ;
    //}
    //if (strcmp(argv[1], "bedsearch") == 0) {
    //    load_chromosomes(c->reference) ;
    //    std::ifstream bed_file(c->bed) ;
    //    std::string line ;
    //    vector<string> lines ;
    //    while (std::getline(bed_file, line)) {
    //        istringstream iss(line) ;
    //        vector<string> tokens{istream_iterator<string>{iss}, istream_iterator<string>{}} ;
    //        lines.push_back(tokens[3]) ;
    //    }
    //    vector<int> n ;
    //    vector<int> m ;
    //    for(int i = 0; i < 48; i++) {
    //        n.push_back(0) ;
    //        m.push_back(0) ;
    //    }
    //    cout << "Searching for " << lines.size() << " sequences.." << endl ;
    //    #pragma omp parallel for num_threads(48)
    //    for (int i = 0; i < lines.size(); i++) {
    //        int t = omp_get_thread_num() ;
    //        char* s = nullptr ;
    //        for (auto chrom = chromosome_seqs.begin(); chrom != chromosome_seqs.end(); chrom++) {
    //            s = strstr(chromosome_seqs[chrom->first], lines[i].c_str()) ;
    //            if (s != nullptr) {
    //                m[t] += 1 ;
    //                cout << "Match " << lines[i] << endl ;
    //                break ;
    //            }
    //        }
    //        n[t] += 1 ;
    //    }
    //    int _n ;
    //    int _m ;
    //    for(int i = 0; i < 48; i++) {
    //        _n += n[i] ;
    //        _m += m[i] ;
    //    }
    //    cout << "Found " << _m << " out of " << _n << endl ;
    //    exit(0) ;
    //}
    if (argc == 1) {
        cerr << "Usage: " << endl;
        cerr << "\tTo index sample:" << endl ;
        cerr << "stella pingpong index [--binary] [--append /path/to/binary/index] --fastq /path/to/fastq/file [--threads threads] --workdir /outpout/directory" << endl ;
        cerr << "\t\tOptional arguments: " << endl ;
        cerr << "\t\t\t-b, --binary          output index in binary format" << endl ;
        cerr << "\t\t\t-a, --append          append to existing index (must be stored in binary)" << endl ;
        cerr << "\t\t\t-t, --threads         number of threads (default is 1)" << endl ;
        cerr << "\tTo search for specific strings:" << endl ;
        cerr << "stella pingpong search [--index /path/to/index] [--fastq /path/to/fastq] [--threads threads] --workdir /output/directory" << endl ;
        cerr << "\t\tOptional arguments: " << endl ;
        cerr << "\t\t\t--aggregate         aggregate ouputs directly." << endl ;
        cerr << "\t\t\t--cutof             sets cutoff for minimum string abundance (tau)" << endl ;
        cerr << "\tTo aggregate specfici strings:" << endl ;
        cerr << "\t\tstella aggregate --workdir /path/to/string/batches --threads <threads> --cutoff <minimum abundance for strings>" << endl ;
        exit(0) ;
    }
    if (strcmp(argv[1], "scan") == 0) {
        c->parse(argc - 1, argv + 1) ;
        create_workdir() ;
        auto scanner = new Scanner() ;
        scanner->run() ;
    }
    if (strcmp(argv[1], "find") == 0) {
        c->parse(argc - 1, argv + 1) ;
        create_workdir() ;
        auto finder = new Finder() ;
        finder->run() ;
    }
    if (strcmp(argv[1], "extract") == 0) {
        c->parse(argc - 1, argv + 1) ;
        create_workdir() ;
        auto extractor = new Extractor() ;
        extractor->run() ;
    }
    if (strcmp(argv[1], "pingpong") == 0) {
        c->parse(argc - 2, argv + 2) ;
        create_workdir() ;
        auto pingpong = new PingPong() ;
        if (strcmp(argv[2], "index") == 0) {
            pingpong->index() ;
        }
        if (strcmp(argv[2], "search") == 0) {
            pingpong->search() ;
            if (c->aggregate) {
                auto aggregator = new Aggregator() ;
                aggregator->num_batches = pingpong->num_output_batches - 1 ;
                aggregator->run() ;
            }
        }
        exit(0) ;
    }
    if (strcmp(argv[1], "aggregate") == 0) {
        c->parse(argc - 1, argv + 1) ;
        create_workdir() ;
        auto aggregator = new Aggregator() ;
        aggregator->run() ;
    }
}

