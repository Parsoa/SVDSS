#ifndef SV_HPP
#define SV_HPP

#include <string>
#include <iostream>

using namespace std;

class SV
{
public:
    string type;
    string chrom;
    int s;
    string refall;
    string altall;
    int e;
    int l;
    int w;
    int cov;
    int nv;
    int score;
    string gt;
    bool imprecise;

    SV();
    SV(const string type_, const string &chrom_, uint s_, const string &refall_, const string &altall_, const uint w_, const uint cov_, int l_, int nv_, int score_, bool imprecise_ = false);
    void genotype();

    friend ostream &operator<<(ostream &os, const SV &sv);
};

#endif