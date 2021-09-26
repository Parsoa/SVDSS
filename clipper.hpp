#ifndef CLIPPER_HPP
#define CLIPPER_HPP

#include <map>
#include <set>
#include <vector>
#include <string>

#include "htslib/sam.h"
#include "htslib/hts.h"
#include "interval_tree.hpp"

#include "sv.hpp"

using namespace lib_interval_tree;

struct Clip {
    std::string name ;
    uint p ;
    uint l ;
    bool starting ;
    uint w ;

    Clip() { 
        w = 0;
    }

    Clip(const std::string &name_, uint p_, uint l_, bool starting_, uint w_ = 0) {
        name = name_;
        p = p_;
        l = l_;
        starting = starting_;
        w = w_;
    }

    bool operator<(const Clip &c) const {
        return p < c.p;
    }
};

class Clipper {

private:
    std::string chrom ;

    samFile *sfs_bam ;
    bam_hdr_t *sfs_bamhdr ;
    hts_idx_t *sfs_bamindex ;

    std::vector<Clip> clips ;
    std::vector<Clip> extract_clips();
    std::vector<Clip> remove_duplicates(const std::vector<Clip> &);
    std::vector<Clip> combine(const std::vector<Clip> &);
    std::vector<Clip> filter_lowcovered(const std::vector<Clip> &, const uint);
    std::vector<Clip> cluster(const std::vector<Clip> &, uint);
    std::vector<Clip> filter_tooclose_clips(const std::vector<Clip> &, interval_tree_t<int> &);

public:
    std::vector<SV> svs;

    Clipper(const std::string&, const std::vector<Clip>&) ;
    void call(const std::string &, interval_tree_t<int> &) ;
};

#endif
