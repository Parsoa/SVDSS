#ifndef CALLER_HPP
#define CALLER_HPP

#include <ctime>
#include <dirent.h>
#include <iostream>
#include <map>

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <spdlog/spdlog.h>

#include "chromosomes.hpp"
#include "clipper.hpp"
#include "config.hpp"
#include "extender.hpp"
#include "sfs.hpp"
#include "sv.hpp"
#include "vcf.hpp"

using namespace std;

class Caller {

public:
  void run();

private:
  Configuration *config;
  unordered_map<string, vector<SFS>> SFSs;

  void print_vcf_header();
};

#endif
