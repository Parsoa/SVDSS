#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <iostream>
#include <string>
#include <vector>

#include "bed_utils.hpp"
#include "cxxopts.hpp"

using namespace std;

class Configuration {

private:
  static Configuration *instance;

public:
  static Configuration *getInstance();

  void parse(int argc, char *argv[]);

  int cutoff = 0;
  int overlap = -1;
  int threads = 4;
  int coverage = 50;
  int batch_size = 1000;
  int min_sv_length = 25;
  int min_indel_length = 20;
  int aggregate_batches = 5;
  int min_cluster_weight = 2;
  float min_ratio = 0.97; // FIXME: change name
  float al_accuracy = 0.02;

  bool binary = false;
  bool clipped = false;
  bool putative = true;
  bool assemble = false;
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
  std::string fasta;
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
