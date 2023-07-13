#ifndef CHR_HPP
#define CHR_HPP

#include <fstream>
#include <iomanip>
#include <string>
#include <unordered_map>
#include <vector>

#include <zlib.h>

#include "config.hpp"
#include "kseq.h"

using namespace std;

extern std::vector<std::string> chromosomes;
extern std::unordered_map<std::string, char *> chromosome_seqs;

int get_reference_size(std::ifstream &fasta_file);
void load_chromosomes(std::string path);
std::string canonicalize(std::string);
std::string reverse_complement(std::string);
std::string load_chromosome(string path);

#endif
