#include <assert.h>
#include <cstdint>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <unordered_map>

#include "chromosomes.hpp"
#include "config.hpp"
#include "lprint.hpp"
#include "ping_pong.hpp"

#include "assembler.hpp"
#include "caller.hpp"
#include "smoother.hpp"

#ifdef LOCAL_BUILD
#include "extractor.hpp"
#include "finder.hpp"
#include "haplotype_shifter.hpp"
#include "kmer_finder.hpp"
#include "scanner.hpp"
#include "shifter.hpp"
#endif

using namespace std;

const string version = "v1.0.4";

// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\

void create_workdir() {
  auto c = Configuration::getInstance();
  struct stat info;
  lprint({c->workdir}, 0);
  if (stat(c->workdir.c_str(), &info) != 0) {
    lprint({"Working directory does not exist. creating.."}, 0);
    int check = mkdir(c->workdir.c_str(), 0777);
    if (check != 0) {
      lprint({"Failed to create output directory, aborting.."}, 2);
      exit(check);
    }
  }
}

void print_help() {
  cerr << endl;
  cerr << "Usage: " << endl;
  cerr << "\t* Index reference/sample:" << endl;
  cerr << "\t\tSVDSS index --fastq/--fasta /path/to/genome/file --index "
          "/path/to/output/index/file"
       << endl;
  cerr << endl;
  cerr << "\t\tOptional arguments: " << endl;
  cerr << "\t\t\t-b, --binary\t\t\t\toutput index in binary format. Allows for "
          "another index to be appended to this index later."
       << endl;
  cerr << "\t\t\t-a, --append /path/to/binary/index\tappend to existing binary "
          "index."
       << endl;
  cerr << endl;
  cerr << "\t* Extract SFS from BAM/FASTQ/FASTA files:" << endl;
  cerr << "\t\tSVDSS search --index /path/to/index --fastq/--bam "
          "/path/to/input --workdir /output/directory"
       << endl;
  cerr << endl;
  cerr << "\t\tOptional arguments: " << endl;
  cerr << "\t\t\t--assemble\t\t\t\tautomatically runs SVDSS assemble on output"
       << endl;
  cerr << endl;
  cerr << "\t* Assemble SFS into superstrings:" << endl;
  cerr << "\t\tSVDSS assemble --workdir /path/to/.sfs/files --batches "
          "/number/of/SFS/batches"
       << endl;
  cerr << endl;
  cerr << "\t* Reconstruct sample:" << endl;
  cerr << "\t\tSVDSS smooth --workdir /output/file/direcotry --bam "
          "/path/to/input/bam/file --reference /path/to/reference/genome/fasta"
       << endl;
  cerr << endl;
  cerr << "\t* Call SVs:" << endl;
  cerr << "\t\tSVDSS call --workdir /path/to/assembled/.sfs/files --bam "
          "/path/to/input/bam/file --reference /path/to/reference/genome/fasta"
       << endl;
  cerr << endl;
  cerr << "\t\tOptional arguments: " << endl;
  cerr << "\t\t\t--clipped\t\t\t\tcalls SVs from clipped SFS." << endl;
  cerr << "\t\t\t--min-cluster-weight\t\t\tminimum number of supporting "
          "superstrings for a call to be reported."
       << endl;
  cerr << "\t\t\t--min-sv-length\t\t\t\tminimum length of reported SVs. "
          "Default is 25. Values < 25 are ignored."
       << endl;
  cerr << endl;
  cerr << "General options: " << endl;
  cerr << "\t--threads\t\t\t\t\t\tsets number of threads, default 4." << endl;
  cerr << "\t--version\t\t\t\t\t\tprint version information." << endl;
  cerr << "\t--help\t\t\t\t\t\t\tprint this help message." << endl;
  cerr << endl;
}

int main(int argc, char **argv) {
  cerr << "SVDSS, Structural Variant Discovery from Sample-specific Strings."
       << endl;
  time_t t;
  time(&t);
  auto c = Configuration::getInstance();

  if (argc == 1) {
    print_help();
    exit(1);
  }

  c->parse(argc, argv);

  if (c->version) {
    cout << "SVDSS, " << version << endl;
    exit(0);
  }

  if (argc == 1 || c->help) {
    print_help();
    exit(0);
  }
  cerr << "Mode: " << argv[1] << endl;
  create_workdir();
  if (strcmp(argv[1], "call") == 0) {
    auto caller = new Caller();
    caller->run();
  } else if (strcmp(argv[1], "index") == 0) {
    auto pingpong = new PingPong();
    pingpong->index();
  } else if (strcmp(argv[1], "search") == 0) {
    auto pingpong = new PingPong();
    pingpong->search();
  } else if (strcmp(argv[1], "assemble") == 0) {
    auto assembler = new Assembler();
    assembler->run();
  } else if (strcmp(argv[1], "smooth") == 0) {
    auto smoother = new Smoother();
    smoother->run();
  } else {
    print_help();
    exit(1);
  }
  time_t s;
  time(&s);
  lprint({"Complete. Runtime:", to_string(s - t) + " seconds."});

  return 0;
}
