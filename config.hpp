#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <iostream>
#include <string>
#include <vector>

#include "cxxopts.hpp"

using namespace std;

// clang-format off

static const char MAIN_USAGE_MESSAGE[] =
  "SVDSS [index|smooth|search|call] --help\n"
  "      --help                           print help message\n"
  "      --version                        print version\n";

static const char INDEX_USAGE_MESSAGE[] =
  "SVDSS index --reference <reference> --index <index>\n"
  "      --binary                         store the index in binary format\n"
  "      --append <oldindex>              append to existing index\n"
  "      --help                           print help message\n";

static const char SMOOTH_USAGE_MESSAGE[] =
  "SVDSS smooth --reference <reference> --bam <bam>\n"
  "      --min-mapq                       minimum mapping quality (default: 20)\n"
  "      --threads <INT>                  number of threads to use (default: 4)\n"
  "      --help                           print help message\n";

static const char SEARCH_USAGE_MESSAGE[] =
  "SVDSS search --index <index> [--bam <bam>|--fastq <fastq>]\n"
  "      --bsize <INT>                    batch size (default: 10000)\n"
  "      --noputative                     when input is smoothed bam, do not filter unsmoothed reads (default: putative)\n"
  "      --noassemble                     do not assemble specific strings overlapping on a read (default: assemble)\n"
  "      --threads <INT>                  number of threads to use (default: 4)\n"
  "      --help                           print help message\n";

static const char CALL_USAGE_MESSAGE[] =
  "SVDSS call --reference <reference> --bam <bam> --sfs <sfs>\n"
  "      --poa <FILE>                     store POA in .sam format to this file (default: do not store)\n"
  "      --clusters <FILE>                store clusters to this file (default: do not store)\n"
  "      --min-cluster-weight <INT>       minimum number of supporting superstrings for a call to be reported (default: 2)\n"
  "      --min-sv-length <INT>            minimum length of reported SVs (default: 25)\n"
  "      --noht                           do not use haplotagging information even if present\n"
  "      --noref                          do not report 0/0 calls\n"
  "      --min-mapq                       minimum mapping quality (default: 20)\n"
  "      --clipped                        calls SVs from clipped SFS (EXPERIMENTAL)\n"
  "      --threads <INT>                  number of threads to use (default: 4)\n"
  "      --help                           print help message\n";

// clang-format on

class Configuration {

private:
  static Configuration *instance;

public:
  static Configuration *getInstance();

  void parse(int argc, char *argv[]);
  void print_help(const string &) const;

  // general
  int threads = 4;
  int batch_size = 10000;
  bool version = false;
  bool verbose = false;
  bool help = false;

  // pingpong.index
  bool binary = false;
  // pingpong.search
  bool assemble = true;
  bool putative = true;
  int overlap = -1;
  int max_output = 100000;
  // call
  uint flank = 100;
  uint ksize = 7;
  uint min_sv_length = 25;
  uint min_mapq = 20;
  int min_indel_length = 20;
  uint min_cluster_weight = 2;
  float min_ratio = 0.97; // FIXME: change name
  bool useht = true;
  bool noref = false;
  bool clipped = false;

  string bam = "";
  string sfs = "";
  string append = "";
  string index = "";
  string reference = "";
  string fastq = "";
  string poa = "";
  string clusters = "";

private:
  Configuration();

  Configuration(Configuration const &) = delete;
  void operator=(Configuration const &) = delete;

  Configuration &operator[](string);

  cxxopts::Options parser;
};

#endif
