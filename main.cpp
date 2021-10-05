#include <omp.h>
#include <ctime>
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
#include "lprint.hpp"
#include "ping_pong.hpp"
#include "chromosomes.hpp"

#include "tau.hpp"
#include "caller.hpp"
#include "assembler.hpp"
#include "converter.hpp"
#include "aggregator.hpp"
#include "reconstructor.hpp"

#ifdef LOCAL_BUILD
#include "finder.hpp"
#include "shifter.hpp"
#include "scanner.hpp"
#include "extractor.hpp"
#include "kmer_finder.hpp"
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
    lprint({c->workdir}, 0);
    if (stat(c->workdir.c_str(), &info) != 0) {
      lprint({"Working directory does not exist. creating.."}, 0);
        int check = mkdir(c->workdir.c_str(), 0777) ;
        if (check != 0) {
            lprint({"Failed to create output directory, aborting.."}, 2);
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
    cerr << "Ping-pong, comparative genome analysis using sample-specific string detection in accurate long reads." << endl ;
    time_t t ;
    time(&t) ;
    auto c = Configuration::getInstance() ;
    if (argc == 1) {
        print_help() ;
        exit(0) ;
    }
    c->parse(argc - 1, argv + 1) ;
    create_workdir() ;
    //cerr << "Running on " << c->threads << " threads." << endl ;
    #ifdef LOCAL_BUILD
    if (strcmp(argv[1], "find") == 0) {
        auto finder = new Finder() ;
        finder->run() ;
        exit(0) ;
    }
    if (strcmp(argv[1], "kmer") == 0) {
        auto finder = new KmerFinder() ;
        finder->run() ;
        exit(0) ;
    }
    if (strcmp(argv[1], "extract") == 0) {
        auto extractor = new Extractor() ;
        extractor->run() ;
        exit(0) ;
    }
    if (strcmp(argv[1], "shift-bed") == 0) {
        auto shifter = new HaplotypeShifter() ;
        shifter->shift_bed_file() ;
        exit(0) ;
    }
    if (strcmp(argv[1], "haplotype-shifter") == 0) {
        auto shifter = new HaplotypeShifter() ;
        shifter->load_tracks() ;
        exit(0) ;
    }
    #endif
    if (strcmp(argv[1], "tau") == 0) {
        auto tau = new Tau() ;
        tau->run() ;
        exit(0) ;
    } else if (strcmp(argv[1], "call") == 0) {
        auto caller = new Caller() ;
        caller->run() ;
    } else if (strcmp(argv[1], "index") == 0) {
        auto pingpong = new PingPong() ;
        pingpong->index() ;
    } else if (strcmp(argv[1], "query") == 0) {
        auto pingpong = new PingPong() ;
        bool b = pingpong->query(string(argv[2])) ;
        cerr << (b ? "SFS" : "Not SFS") << endl ;
    } else if (strcmp(argv[1], "search") == 0) {
        auto pingpong = new PingPong() ;
        pingpong->search() ;
        if (c->aggregate) {
            auto aggregator = new Aggregator() ;
            c->aggregate_batches = pingpong->num_output_batches ;
            aggregator->run() ;
        }
    } else if (strcmp(argv[1], "convert") == 0) {
        auto converter = new Converter() ;
        converter->run() ;
    } else if (strcmp(argv[1], "assemble") == 0) {
        auto assembler = new Assembler() ;
        assembler->run() ;
    } else if (strcmp(argv[1], "aggregate") == 0) {
        auto aggregator = new Aggregator() ;
        aggregator->run() ;
    } else if (strcmp(argv[1], "reconstruct") == 0) {
        auto reconstructor = new Reconstructor() ;
        reconstructor->run() ;
    } else {
        print_help() ;
    }
    time_t s ;
    time(&s) ;
    lprint({"Complete. Runtime:", to_string(s - t) + " seconds."}) ;
}

