#ifndef CLUSTER_HPP
#define CLUSTER_HPP

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>

#include "abpoa.h"

struct Cluster {
    std::string chrom ;
    int s ;
    int e ;
    int cov ;
    std::vector<std::string> seqs ;

    Cluster(const std::string &chrom_, uint s_, uint e_, uint cov_ = 0) {
        chrom = chrom_ ;
        s = s_ ;
        e = e_ ;
        cov = cov_ ;
    }

    void set_cov(uint cov_) {
        cov = cov_ ;
    }

    void add(const std::string &seq) {
        seqs.push_back(seq);
    }

    int get_len() const {
        uint l = 0;
        uint n = 0;
        for (const std::string &seq : seqs) {
            ++n;
            l += seq.size();
        }
        return l / n;
    }

    std::vector<std::string> get_seqs() const {
        return seqs ;
    }

    uint size() const {
        return seqs.size() ;
    }

    void dump(std::ofstream &o) const ;
    
    std::string poa() ;
};

#endif
