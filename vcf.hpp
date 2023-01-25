#ifndef VCF_HPP
#define VCF_HPP

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <math.h>
#include <sstream>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <unordered_map>
#include <vector>

#include "lprint.hpp"

using namespace std;

// CHECKME: I created the SV class. Maybe we can merge it with this struct
struct vcf_variant_t {
  std::string chrom;
  int pos;
  std::string ref;
  std::string alleles[2];
  int svlen;

  bool operator==(const vcf_variant_t &v) const {
    return chrom == v.chrom && pos == v.pos;
  }
};

namespace std {
template <> struct hash<vcf_variant_t> {
  std::size_t operator()(const vcf_variant_t &v) const {
    return std::hash<std::string>()(v.chrom) + std::hash<int>()(v.pos);
  }
};
} // namespace std

std::unordered_map<std::string, std::vector<vcf_variant_t>>
    load_vcf_file(std::string);

#endif
