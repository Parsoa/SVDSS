#ifndef EXTENDER_HPP
#define EXTENDER_HPP

#include <set>
#include <omp.h>
#include <string>
#include <vector>
#include <fstream>
#include <stdlib.h>

#include "parasail.h"
#include "htslib/sam.h"
#include "htslib/hts.h"
#include "interval_tree.hpp"
#include "parasail/matrices/blosum62.h"

#include "sv.hpp"
#include "bam.hpp"
#include "sfs.hpp"
#include "config.hpp"
#include "cluster.hpp"
#include "clipper.hpp"
#include "chromosomes.hpp"

using namespace lib_interval_tree;

class Extender {

private:
    
    uint minw = 2 ;
    uint mind = 15 ;
    uint skip_1 = 0 ;      // SFS skipped since no first/last base can be placed from read alignment (should be rare)
    uint skip_2 = 0 ;      // SFS skipped since it couldn't be extended
    uint skip_3 = 0 ;      // SFS skipped since reads starts/ends inside a cluster
    uint small_cl = 0 ;    // number of cluster (before clustering) with low support
    uint extcl = 0 ;       // number of extended clusters (after clustering)
    uint small_extcl = 0 ; // number of extended clusters (after clustering) with low support
    uint maxw = 100 ;
    uint minsupp = 2 ;
    uint kmer_size = 7 ;

    Configuration* config ;

    samFile* bam_file ;
    hts_idx_t* bam_index ;
    bam_hdr_t* bam_header ;

    std::vector<Cluster> clusters ;
    std::vector<ExtCluster> ext_clusters ;
    std::unordered_map<std::string, std::vector<SFS>>* SFSs ;
    std::vector<ExtSFS> extended_sfs ;
    
    void extend_parallel() ;
    void extend_alignment(bam1_t* aln, int index) ;
    void process_batch(vector<bam1_t*> bam_entries, int) ;
    bool load_batch_bam(int threads, int batch_size, int p) ;
    std::pair<int, int> get_unique_kmers(const std::vector<std::pair<int, int>> &alpairs, const uint k, const bool from_end, std::string chrom) ;

    void cluster() ;
    void cluster_interval_tree() ;

    void call() ;
    std::vector<std::pair<uint, char>> parse_cigar(std::string) ;
    std::vector<Cluster> cluster_by_length(const Cluster& cluster) ;
    
    // parallelize
    int threads ; 
    int batch_size ;
    std::vector<std::vector<SV>> _p_svs ;
    std::vector<std::vector<Clip>> _p_clips ;
    std::vector<std::vector<ExtSFS>> _p_extended_sfs ;
    std::vector<std::vector<Cluster>> _p_clusters ;
    std::vector<std::vector<Consensus>> _p_alignments ;
    std::vector<std::vector<std::vector<bam1_t*>>> bam_entries ;
    std::vector<std::unordered_map<string, interval_tree_t<int>>> _p_tree ;
    std::vector<std::map<std::pair<int, int>, std::vector<ExtSFS>>> _p_sfs_clusters ;

public:
    Extender(std::unordered_map<std::string, std::vector<SFS>>*) ;
    
    std::vector<SV> svs ;
    std::vector<Clip> clips ;
    std::vector<Consensus> alignments ;

    void run(int threads) ;
};

#endif
