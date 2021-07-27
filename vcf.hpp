#ifndef VCF_HPP
#define VCF_HPP

#include <vector>
#include <string>
#include <iterator>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <math.h>
#include <string.h>

#include "lprint.hpp"

using namespace std ;

// CHECKME: I created the SV class. Maybe we can merge it with this struct
struct vcf_variant_t {
    std::string chrom ;
    int pos ;
    std::string ref ;
    std::string alleles[2] ;
    int svlen ;
    
    bool operator==(const vcf_variant_t& v) const {
        return chrom == v.chrom && pos == v.pos ; 
    }
};

namespace std {
    template <> struct hash<vcf_variant_t> {
        std::size_t operator()(const vcf_variant_t& v) const {
            return std::hash<std::string>()(v.chrom) + std::hash<int>()(v.pos) ;
        }
    };
}

std::unordered_map<std::string, std::vector<vcf_variant_t>> load_vcf_file(std::string);

// Write custom VCF header to ostream (it needs chromosome_seqs)
void print_vcf_header(const unordered_map<string, char*> &, ofstream &, const string &sample = "DEFAULT");

#endif
