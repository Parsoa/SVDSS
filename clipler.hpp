#ifndef CLIPLER_HPP
#define CLIPLER_HPP

#include <string>
#include <list>
#include <map>
#include <set>

#include "htslib/sam.h"
#include "htslib/hts.h"
#include "interval_tree.hpp"

#include "sv.hpp"

using namespace std;
using namespace lib_interval_tree;

struct Clip
{
    string name;
    uint p;
    uint l;
    bool starting;
    uint w;

    Clip() { w = 0; }

    Clip(const string &name_, uint p_, uint l_, bool starting_, uint w_ = 0)
    {
        name = name_;
        p = p_;
        l = l_;
        starting = starting_;
        w = w_;
    }

    bool operator<(const Clip &c) const
    {
        return p < c.p;
    }
};

class Clipler {

private:
    string chrom;
    samFile *sfs_bam;
    bam_hdr_t *sfs_bamhdr;
    hts_idx_t *sfs_bamindex;

    list<Clip> extract_clips();
    list<Clip> remove_duplicates(const list<Clip> &);
    list<Clip> combine(const list<Clip> &);
    list<Clip> filter_lowcovered(const list<Clip> &, const uint);
    list<Clip> cluster(const list<Clip> &, uint);
    list<Clip> filter_tooclose_clips(const list<Clip> &, interval_tree_t<int> &);

public:
    list<SV> osvs;

    Clipler(const string&, samFile *, bam_hdr_t *, hts_idx_t *);
    void call(const string &, interval_tree_t<int> &);
};

#endif
