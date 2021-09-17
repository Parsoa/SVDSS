#ifndef INSDELLER_HPP
#define INSDELLER_HPP

#include <set>
#include <omp.h>
#include <string>
#include <vector>
#include <fstream>
#include <stdlib.h>

#include "ksw2.h"
#include "edlib.hpp"
#include "htslib/sam.h"
#include "htslib/hts.h"
#include "interval_tree.hpp"
#include "rapidfuzz/fuzz.hpp"
#include "rapidfuzz/utils.hpp"

#include "sv.hpp"
#include "bam.hpp"
#include "config.hpp"
#include "cluster.hpp"
#include "haplotyper.hpp"

using namespace lib_interval_tree;

class Insdeller {

private:
    std::string chrom;

    samFile* sfs_bam ;
    bam_hdr_t* sfs_bamhdr ;
    hts_idx_t* sfs_bamindex ;

    std::vector<samFile*> read_bam;
    std::vector<bam_hdr_t*> read_bamhdr;
    std::vector<hts_idx_t*> read_bamindex;

    double expected_mismatch_rate ;

    Haplotyper haplotyper ;
    Configuration* config ;

    // Cluster superstrings by mapping position
    std::vector<Cluster> position_cluster();
    Cluster merge_close_fragments(const Cluster& cluster, int distance) ;
    // Type clustering - 2 clusters with all fragments with only I and only D or 1 cluster with both (if we have at least one mixed fragment)
    std::vector<Cluster> type_cluster(const Cluster &);
    // Fragments are extended + Fragments on same reads are merged
    Cluster extend(Cluster &, int);
    // SequenceSimilarity-based clustering
    std::vector<Cluster> scluster(const Cluster &);
    // Extract SVs (poa + realign)
    Cluster compress_cluster(const Cluster& c, int size) ;
    std::vector<Cluster> cluster_breakpoints(const Cluster& cluster, float ratio) ;
    std::vector<SV> call_batch(std::vector<Cluster>& position_clusters, const string&, ofstream&) ;
    std::vector<SV> call_svs(const Cluster& cluster, const string&) ;
    std::vector<SV> call_poa_svs(Cluster&, const string&, std::ofstream &o);
    std::vector<SV> filter_chain_svs(std::vector<SV> svs) ;
    bool should_filter_read(bam1_t* alignment, char* read_seq, string chrom, int* global_num_bases, int* global_num_mismatch) ;
    // Global realignment of consensus and subreference
    CIGAR align_ksw2(const char *, const char *, int, int, int, int);
    CIGAR align_edlib(const char *, const char *, int, int, int, int);
    // Remove svs that are duplicate (same variant on different representative)
    std::vector<SV> remove_duplicate_svs(const vector<SV> &);
    
    int cluster_anchors(const Cluster&) ;

public:
    std::vector<SV> osvs;
    interval_tree_t<int> vartree; // this will contain the breakpoints of all (precise) variations that are considered for output

    Insdeller(const std::string&);
    void call(const std::string &, ofstream &);
};

#endif
