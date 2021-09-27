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

    std::string chrom ;
    std::string ref_seq ;
    Configuration* config ;

    samFile* bam_file ;
    bam_hdr_t* bam_header ;
    hts_idx_t* bam_index ;
    hts_itr_t* bam_iter ;

    std::vector<ExtSFS> extended_sfs ;
    std::vector<Cluster> clusters ;
    std::unordered_map<std::string, std::vector<SFS>>* SFSs ;

    void extend_alignment(bam1_t* aln) ;
    
    std::pair<int, int> get_unique_kmers(const std::vector<std::pair<int, int>> &alpairs, const uint k, const bool from_end) ;
    std::vector<Cluster> cluster_by_length(const Cluster& cluster) ;
    std::vector<std::pair<uint, char>> parse_cigar(std::string) ;


    // parallelize
    int threads ; 
    int batch_size ;
    std::vector<std::vector<Clip>> _p_clips ;
    std::vector<std::vector<Consensus>> _p_alignments ;
    std::vector<std::vector<std::vector<bam1_t*>>> bam_entries ;

    void extended_parallel() ;
    void process_batch(vector<bam1_t*> bam_entries) ;
    bool load_batch_bam(int threads, int batch_size, int p) ;

public:
    Extender(const std::string&, std::unordered_map<std::string, std::vector<SFS>>*) ;
    
    std::vector<SV> svs ;
    std::vector<Clip> clips ;
    std::vector<Consensus> alignments ;

    void call() ;
    void extend() ;
    void cluster() ;
};

#endif
