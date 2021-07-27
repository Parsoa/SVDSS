#ifndef INSDELLER_HPP
#define INSDELLER_HPP

#include <fstream>
#include <string>
#include <list>
#include <set>

#include "htslib/sam.h"
#include "htslib/hts.h"
#include "rapidfuzz/fuzz.hpp"
#include "rapidfuzz/utils.hpp"
#include "ksw2.h"
#include "interval_tree.hpp"


#include "bam.hpp"
#include "cluster.hpp"
#include "sv.hpp"

using namespace std;
using namespace lib_interval_tree;

class Insdeller {

private:
    string chrom;

    samFile *sfs_bam;
    bam_hdr_t *sfs_bamhdr;
    hts_idx_t *sfs_bamindex;

    samFile *read_bam;
    bam_hdr_t *read_bamhdr;
    hts_idx_t *read_bamindex;

    // Cluster superstrings by mapping position
    list<Cluster> pcluster();
    // Type clustering - 2 clusters with all fragments with only I and only D or 1 cluster with both (if we have at least one mixed fragment)
    list<Cluster> tcluster(const Cluster &);
    // Fragments are extended + Fragments on same reads are merged
    Cluster extend(const Cluster &);
    // SequenceSimilarity-based clustering
    list<Cluster> scluster(const Cluster &);
    // Extract SVs (poa + realign)
    list<SV> extract(const Cluster &, const string &, ofstream &);
    // Global realignment of consensus and subreference
    CIGAR align(const char *, const char *, int, int, int, int);
    // Merge svs that are on the same representative
    list<SV> merge_svs(const list<SV> &, const string &);
    // Remove svs that are duplicate (same variant on different representative)
    list<SV> dedup_svs(const list<SV> &);

public:
    list<SV> osvs;
    interval_tree_t<int> vartree; // this will contain the breakpoints of all (precise) variations that are considered for output

    Insdeller(const string&, samFile *, bam_hdr_t *, hts_idx_t *, samFile *, bam_hdr_t *, hts_idx_t *);
    void call(const string &, ofstream &);
};

#endif
