#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <iostream>
#include <string>
#include <vector>

#include "bed_utils.hpp"
#include "cxxopts.hpp"

using namespace std;

static const char *MAIN_USAGE_MESSAGE =
    "SVDSS [index|smooth|search|call] --help\n"
    "      -h,--help                          print help message\n"
    "      --version                          print version\n";

static const char *INDEX_USAGE_MESSAGE =
    "SVDSS index [-h] --reference <reference> --index <index>\n"
    "      -b,--binary                        store the index in binary "
    "format\n"
    "      -a,--append <oldindex>             append to existing index\n";

static const char *SMOOTH_USAGE_MESSAGE =
    "SVDSS smooth [-h] --reference <reference> --bam <bam>\n"
    "      --threads                          number of threads to use "
    "(default: 4)\n";
;

static const char *SEARCH_USAGE_MESSAGE =
    "SVDSS search [-h] --index <index> --bam <bam> --workdir <wd>\n"
    "      --threads                          number of threads to use "
    "(default: 4)\n"
    "      --noputative                       when input is smoothed bam, do "
    "not filter unsmoothed reads (default: putative)\n"
    "      --noassemble                       do not assemble specific strings "
    "overlapping on a read (default: assemble)\n";

static const char *CALL_USAGE_MESSAGE =
    "SVDSS call [-h] --reference <reference> --bam <bam> --workdir <wd>\n"
    "      --min-cluster-weight               minimum number of supporting "
    "superstrings for a call to be reported (default: 2)\n"
    "      --min-sv-length                    minimum length of reported SVs "
    "(default: 25)\n"
    "      --threads                          number of threads to use "
    "(default: 4)\n"
    "      --clipped                          calls SVs from clipped SFS "
    "(EXPERIMENTAL)\n";

class Configuration {

private:
  static Configuration *instance;

public:
  static Configuration *getInstance();

  void parse(int argc, char *argv[]);
  void print_help(const string &) const;

  int cutoff = 0;
  int overlap = -1;
  int threads = 4;
  int coverage = 50;
  int batch_size = 1000;
  int min_sv_length = 25;
  int min_indel_length = 20;
  int min_cluster_weight = 2;
  float min_ratio = 0.97; // FIXME: change name
  float al_accuracy = 0.02;

  bool binary = false;
  bool clipped = false;
  bool putative = true;
  bool assemble = true;
  bool aggregate = false;
  bool selective = true;
  bool version = false;
  bool verbose = false;
  bool help = false;

  std::string bed;
  std::string bam;    // reads bam (reconstructed or not)
  std::string sfsbam; // superstrings bam (from realignment)
  std::string vcf;
  std::string type;
  std::string workdir;
  std::string append;
  std::string index;
  std::string fastq;
  std::string target;
  std::string prefix;
  std::string reference;

private:
  Configuration();

  Configuration(Configuration const &) = delete;
  void operator=(Configuration const &) = delete;

  Configuration &operator[](std::string);

  cxxopts::Options parser;
};

#endif
