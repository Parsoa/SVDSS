#include <math.h>

#include "tau.hpp"
#include "chromosomes.hpp"

void Tau::run() {
    config = Configuration::getInstance() ;
    int i = 0 ;
    string inpath = config->workdir + "/superstring_loci.bed" ;
    string line ;
    ifstream in(inpath) ;
    while (getline(in, line)) {
        istringstream iss(line) ;
        vector<string> tokens{istream_iterator<string>{iss}, istream_iterator<string>{}} ;
        auto chr = tokens[2] ;
        auto pos = std::stoi(tokens[3]) ;
        auto sfs = canonicalize(tokens[1]) ;
        tau_distribution[sfs].sfsnames.push_back(tokens[0]) ;
        if (tau_distribution[sfs].loci.find(chr) == tau_distribution[sfs].loci.end()) {
            tau_distribution[sfs].loci[chr].push_back({chr, pos, 0}) ;
            tau_distribution[sfs].num_loci += 1 ;
        } else {
            bool found = false ;
            for (auto& locus: tau_distribution[sfs].loci[chr]) {
                if (abs(locus.position - pos) < 100) {
                    found = true ;
                    break ;
                }
            }
            if (!found) {
                tau_distribution[sfs].loci[chr].push_back({chr, pos, 0}) ;
                tau_distribution[sfs].num_loci += 1 ;
            }
        }
        i++ ;
        if (i % 100000 == 0) {
            cout << "Loaded " << i << " alignments for " << tau_distribution.size() << " SFS." << endl ;
        }
    }
    cout << "Loaded " << i << " alignments for " << tau_distribution.size() << " SFS." << endl ;
    ofstream l_o("tau_distribution.bed") ;
    for (auto &sfs: tau_distribution) {
        int i = 0 ;
        for (auto& sfsname: sfs.second.sfsnames) {
            l_o << (i == 0 ? sfs.first : "*") << "\t" << sfsname << "\t" << to_string(sfs.second.num_loci) << endl ;
            i += 1 ;
        }
    }
}

