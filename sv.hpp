#ifndef SV_HPP
#define SV_HPP

#include <string>
#include <iostream>

using namespace std;

class SV {

public:
    string type ;
    string chrom ;
    string idx ;
    uint s ;
    uint e ;
    string refall ;
    string altall ;
    uint w ;
    uint cov ;
    int l ;
    int ngaps ;
    int score ;
    string gt ;
    bool imprecise ;

    SV();
    SV(const std::string type_, const string &chrom_, uint s_, const string &refall_, const string &altall_, const uint w_, const uint cov_, const int ngaps_, const int score_, bool imprecise_ = false, uint l_ = 0) ;
    void genotype();

    friend ostream &operator<<(ostream &os, const SV &sv) ;
};

#endif
