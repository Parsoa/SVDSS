#include <math.h>
#include <string>
#include <fstream>
#include <sstream>
#include <iterator>
#include <stdlib.h>
#include <iostream>
#include <algorithm>

#include "vcf.hpp"

using namespace std ;

unordered_map<string, vector<vcf_variant_t>> load_vcf_file(string path) {
    std::unordered_map<std::string, std::vector<vcf_variant_t>> vcf_variants ;
    std::ifstream vcf_file(path) ;
    std::string line ;
    unordered_map<string, int> header ;
    int i = 0 ;
    int n = 0 ;
    string chrom = "UNKNOWN" ;
    while (std::getline(vcf_file, line)) {
        if (line[0] == '#' and line[1] == '#') {
            continue ;
        } else if (line[0] == '#') {
            //TODO parse header
            istringstream iss(line) ;
            vector<string> tokens{istream_iterator<string>{iss}, istream_iterator<string>{}} ;
        } else {
            istringstream iss(line) ;
            vector<string> tokens{istream_iterator<string>{iss}, istream_iterator<string>{}} ;
            if (chrom != tokens[0]) {
                chrom = tokens[0] ;
            }
            vcf_variant_t vcf_variant ;
            vcf_variant.chrom = tokens[0] ;
            vcf_variant.pos = std::stoi(tokens[1]) ;
            vcf_variant.ref = tokens[3] ;
            int l = vcf_variants[chrom].size() ;
            if (l != 0 && vcf_variants[chrom][l - 1] == vcf_variant) {
                vcf_variants[chrom][l - 1].alleles[1] = tokens[4] ;
                vcf_variants[chrom][l - 1].svlen = max(vcf_variants[chrom][l - 1].svlen, int(tokens[4].length() - tokens[3].length())) ;
            } else {
                vcf_variant.alleles[0] = tokens[4] ;
                vcf_variant.alleles[1] = "$" ;
                vcf_variants[chrom].push_back(vcf_variant) ;
                vcf_variants[chrom][l - 1].svlen = tokens[4].length() - tokens[3].length() ;
            }
            n++ ;
        }
    }
    cout << "Loaded " << n << " variants from " << path << endl ;
    return vcf_variants ;
}
