#ifndef SV_HPP
#define SV_HPP

#include <map>
#include <string>
#include <iostream>
#include <unordered_map>

class SV {

public:
    std::string type ;
    std::string chrom ;
    std::string idx ;
    int s ;
    int e ;
    std::string refall ;
    std::string altall ;
    uint w ;
    uint cov ;
    int l ;
    int ngaps ;
    int score ;
    std::string gt ;
    bool imprecise ;

    SV();
    SV(const std::string type_, const std::string &chrom_, uint s_, const std::string &refall_, const std::string &altall_, const uint w_, const uint cov_, const int ngaps_, const int score_, bool imprecise_ = false, uint l_ = 0) ;
    void genotype();
    
    bool operator<(const SV& c) const {
        return s < c.s;
    }

    bool operator==(const SV& c) const {
        return chrom == c.chrom and s == c.s and e == c.e and type == c.type;
    }

    friend std::ostream &operator<<(std::ostream &os, const SV &sv) ;
};

class SVCluster {
public:
    std::string chrom ;
    std::string type ;
    int s ;
    int svlen ;
    std::map<SV, int> svs ;

    SVCluster() ;
    SVCluster(SV& sv) {
        chrom = sv.chrom ;
        type = sv.type ;
        s = sv.s ;
        svlen = std::abs(sv.l) ;
        svs[sv] = 1 ;
    }

    void add_sv(SV sv, const int count = 1) ;

    bool operator<(const SVCluster& c) const {
        return s < c.s;
    }

    bool operator==(const SVCluster& c) const {
        return chrom == c.chrom and s == c.s and type == c.type;
    }
};


#endif
