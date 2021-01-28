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
#include "ping_pong.hpp"
#include "aggregator.hpp"
#include "chromosomes.hpp"

#ifdef LOCAL_BUILD
#include "finder.hpp"
#include "shifter.hpp"
#include "scanner.hpp"
#include "extractor.hpp"
#include "haplotype_shifter.hpp"
#endif

using namespace std ;

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

void print_help() {
    cerr << "Usage: " << endl;
    cerr << "\tTo index sample:" << endl ;
    cerr << "\t\tPingPong index [--binary] [--append /path/to/binary/index] --fastq /path/to/fastq/file [--threads threads] --index /path/to/output/index/file" << endl ;
    cerr << "\t\tOptional arguments: " << endl ;
    cerr << "\t\t\t-b, --binary          output index in binary format" << endl ;
    cerr << "\t\t\t-a, --append          append to existing index (must be stored in binary). DON'T pass this option for building an index you want to use directly." << endl ;
    cerr << "\t\t\t-t, --threads         number of threads (default is 1)" << endl ;
    cerr << "\tTo search for specific strings:" << endl ;
    cerr << "\t\tPingPong search [--index /path/to/index] [--fastq /path/to/fastq] [--threads threads] --workdir /output/directory" << endl ;
    cerr << "\t\tOptional arguments: " << endl ;
    cerr << "\t\t\t--aggregate         aggregate ouputs directly." << endl ;
    cerr << "\t\t\t--cutof             sets cutoff for minimum string abundance (tau)" << endl ;
    cerr << "\tTo aggregate specfic strings:" << endl ;
    cerr << "\t\tPingPong aggregate --workdir /path/to/string/batches --threads <threads> --cutoff <minimum abundance for strings> --batches <number of output batches>" << endl ;
}

int main(int argc, char** argv) {
    cout << "Ping-pong, comparative genome analysis using sample-specific string detection in accurate long reads." << endl ;
    auto c = Configuration::getInstance() ;
    if (argc == 1) {
        print_help() ;
        exit(0) ;
    }
    #ifdef LOCAL_BUILD
    if (strcmp(argv[1], "scan") == 0) {
        c->parse(argc - 1, argv + 1) ;
        create_workdir() ;
        auto scanner = new Scanner() ;
        scanner->run() ;
        exit(0) ;
    }
    if (strcmp(argv[1], "find") == 0) {
        c->parse(argc - 1, argv + 1) ;
        create_workdir() ;
        auto finder = new Finder() ;
        finder->run() ;
        exit(0) ;
    }
    if (strcmp(argv[1], "extract") == 0) {
        c->parse(argc - 1, argv + 1) ;
        create_workdir() ;
        auto extractor = new Extractor() ;
        extractor->run() ;
        exit(0) ;
    }
    if (strcmp(argv[1], "haplotype-shifter") == 0) {
        c->parse(argc - 1, argv + 1) ;
        auto shifter = new HaplotypeShifter() ;
        shifter->run() ;
        exit(0) ;
    }
    #endif
    if (strcmp(argv[1], "index") == 0) {
        c->parse(argc - 1, argv + 1) ;
        create_workdir() ;
        auto pingpong = new PingPong() ;
        pingpong->index() ;
    }
    else if (strcmp(argv[1], "search") == 0) {
        c->parse(argc - 1, argv + 1) ;
        create_workdir() ;
        auto pingpong = new PingPong() ;
        pingpong->search() ;
        if (c->aggregate) {
            auto aggregator = new Aggregator() ;
            c->aggregate_batches = pingpong->num_output_batches ;
            aggregator->run() ;
        }
    }
    else if (strcmp(argv[1], "aggregate") == 0) {
        c->parse(argc - 1, argv + 1) ;
        create_workdir() ;
        auto aggregator = new Aggregator() ;
        aggregator->run() ;
    } else {
        print_help() ;
    }
}

