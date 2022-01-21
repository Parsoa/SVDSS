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
    cerr << "\t\tSVDSS index --fastq/--fasta /path/to/genome/file --index /path/to/output/index/file" << endl ;
    cerr << "\t\tOptional arguments: " << endl ;
    cerr << "\t\t\t-b, --binary\t\t\t\toutput index in binary format. Allows for another index to be appended to this index later." << endl ;
    cerr << "\t\t\t-a, --append /path/to/binary/index\t\t\t\tappend to existing binary index." << endl ;
    cerr << "\tExtract SFS from BAM/FASTQ/FASTA files:" << endl ;
    cerr << "\t\tSVDSS search --index /path/to/index --fastq/--bam /path/to/input --workdir /output/directory" << endl ;
    cerr << "\t\tOptional arguments: " << endl ;
    cerr << "\t\t\t--assemble\t\t\t\tautomatically runs SVDSS assemble on output" << endl ;
    cerr << "\tAssmble SFS into superstrings" << endl ;
    cerr << "\t\tSVDSS assemble --workdir /path/to/.sfs/files --batches /number/of/SFS/batches" << endl ;
    cerr << "\tReconstruct sample" << endl ;
    cerr << "\t\tSVDSS reconstruct --workdir /output/file/direcotry --bam /path/to/input/bam/file --reference /path/to/reference/genome/fasta" << endl ;
    cerr << "\tCall SVs:" << endl ;
    cerr << "\t\tSVDSS call --workdir /path/to/assembled/.sfs/files --bam /path/to/input/bam/file --reference /path/to/reference/genome/fasta" << endl ;
    cerr << "\t\tOptional arguments: " << endl ;
    cerr << "\t\t\t--clipped\t\t\t\tcalls SVs from clipped SFS." << endl ;
    cerr << "\t\t\t--min-cluster-weight\t\t\t\tminimum number of supporting superstrings for a call to be reported." << endl ;
    cerr << "\t\t\t--min-sv-length\t\t\t\tminimum length of reported SVs. Default is 25. Values < 25 are ignored." << endl ;
    cerr << "\tGeneral options: " << endl ;
    cerr << "\t\t--threads\t\t\t\tsets number of threads, default 4." << endl ;
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
    if (strcmp(argv[1], "call") == 0) {
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

