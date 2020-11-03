#ifndef CHR_HPP
#define CHR_HPP

#include <vector>
#include <string>
#include <fstream>
#include <unordered_map>

#include "config.hpp"

extern std::vector<std::string> chromosomes ;
extern std::unordered_map<std::string, char*> chromosome_seqs ;

int get_reference_size(std::ifstream &fasta_file) ;
void load_chromosomes(std::string path) ;
std::string canonicalize(std::string) ;
std::string reverse_complement(std::string) ;

#endif
