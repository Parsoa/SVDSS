#ifndef CALLER_HPP
#define CALLER_HPP

#include <iostream>

#include "htslib/hts.h"
#include "htslib/sam.h"

#include "chromosomes.hpp"
#include "clipper.hpp"
#include "config.hpp"
#include "extender.hpp"
#include "sv.hpp"
#include "vcf.hpp"

class Caller {

public:
  void run();

private:
  Configuration *config;
  std::unordered_map<std::string, std::vector<SFS>> SFSs;

  ofstream ovcf;
  ofstream osam;
  void load_input_sfs();
  void print_vcf_header();
};

#endif
