#ifndef VCF_HPP
#define VCF_HPP

#include <vector>
#include <string> 
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <unordered_map>

struct vcf_variant_t {
    std::string chrom ;
    int pos ;
    std::string ref ;
    std::string alleles[2] ;
    
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

#endif
