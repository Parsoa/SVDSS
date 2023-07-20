#ifndef CHR_HPP
#define CHR_HPP

#include <fstream>
#include <iomanip>
#include <string>
#include <unordered_map>
#include <vector>

#include <zlib.h>

#include <spdlog/spdlog.h>

#include "config.hpp"
#include "kseq.h"

using namespace std;

extern vector<string> chromosomes;
extern unordered_map<string, char *> chromosome_seqs;

void load_chromosomes(string path);

#endif
