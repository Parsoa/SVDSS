#include <map>

#include "insdeller.hpp"

#define ANCHOR_SIZE 7

#define EXTENSION_TYPE_FULL -2
#define EXTENSION_TYPE_KMER -1
#define EXTENSION_TYPE_NONE 0

struct pair_hash {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1,T2> &p) const {
        auto h1 = std::hash<T1>{}(p.first) ;
        auto h2 = std::hash<T2>{}(p.second) ;
        return h1 ^ h2 ;
    }
};

Insdeller::Insdeller(const string &chrom_) {
    chrom = chrom_;
    config = Configuration::getInstance() ;
    sfs_bam = hts_open(config->sfsbam.c_str(), "r") ;
    sfs_bamhdr = sam_hdr_read(sfs_bam) ;
    sfs_bamindex = sam_index_load(sfs_bam, config->sfsbam.c_str()) ;
    expected_mismatch_rate = 0.005 ;
    for (int i = 0; i < config->threads; i++) {
        read_bam.push_back(hts_open(config->bam.c_str(), "r")) ;
        read_bamhdr.push_back(sam_hdr_read(read_bam[i])) ;
        read_bamindex.push_back(sam_index_load(read_bam[i], config->bam.c_str())) ;
    }
}

// cluster SFS based on position
vector<Cluster> Insdeller::position_cluster() {
    lprint({"Positioning reads into clusters.."}) ;
    bam1_t *aln = bam_init1();
    hts_itr_t *itr = sam_itr_querys(sfs_bamindex, sfs_bamhdr, chrom.c_str());

    vector<Cluster> clusters; // this will store the clusters
    clusters.push_back(Cluster(chrom));

    int buffer_len = 10000 ;
    uint8_t* seq = (uint8_t*) malloc(buffer_len) ;
    uint8_t* qual = (uint8_t*) malloc(buffer_len) ;

    int n = 0 ;
    while (sam_itr_next(sfs_bam, itr, aln) > 0) {
        //
        string sfs_name(bam_get_qname(aln));
        int u = sfs_name.rfind(".") ;
        string qname = sfs_name.substr(0, u) ;
        string coords = sfs_name.substr(u + 1, sfs_name.length() - u - 1) ;
        u = coords.find('-') ;
        int read_s = std::stoi(coords.substr(0, u)) ;
        int read_e = std::stoi(coords.substr(u + 1, coords.length() - u - 1)) ;
        //cout << qname << " " << coords << " " << read_s << "-" << read_e << endl ;
        uint ref_s = aln->core.pos ;
        uint ref_e = bam_endpos(aln) ;

        bool has_i = false ;
        bool has_d = false ;
        uint32_t *cigar = bam_get_cigar(aln);
        for (uint i = 0; i < aln->core.n_cigar; i++) {
            uint op = bam_cigar_op(*(cigar + i)) ;
            uint len = cigar[i] >> 4 ;
            // if they don't have long enough indels thennjust ignore them
            if (len > 10) {
                has_i |= op == BAM_CINS ;
                has_d |= op == BAM_CDEL ;
            }
        }

        if (!has_i && !has_d) {
            continue ;
        }
        uint ft = 0 ;
        if (has_i) {
            if (has_d) {
                ft = 2 ;
            } else {
                ft = 1 ;
            }
        }

        auto aligned_pairs = get_aligned_pairs(aln) ;

        uint32_t l = aln->core.l_qseq ; //length of the read
        if (l >= buffer_len) {
            free(seq) ;
            free(qual) ;
            buffer_len = l + 1 ;
            seq = (uint8_t*) malloc(buffer_len) ;
            qual = (uint8_t*) malloc(buffer_len) ;
        }
        uint8_t *s = bam_get_seq(aln) ; //quality string
        uint8_t *q = bam_get_qual(aln) ;
        for (int i = 0; i < l; i++) {
            seq[i] = seq_nt16_str[bam_seqi(s, i)]; //gets nucleotide id and converts them into IUPAC id.
            qual[i] = 'I' ; //q[i] ;
        }
        seq[l] = '\0' ; // null terminate
        qual[l] = '\0' ;

        Cluster &c = clusters.back();
        Fragment f = Fragment(qname, ref_s, ref_e, read_s, read_e, ft, (char*)seq, (char*)qual) ;
        f.aligned_pairs = aligned_pairs ;
        if (c.empty()) {
            c.add_fragment(f) ;
        }
        else {
            if (c.back().ref_e + 500 > ref_s) {
                c.add_fragment(f) ;
            }
            else {
                clusters.push_back(Cluster(chrom)) ;
                clusters.back().add_fragment(f) ;
            }
        }
        n += 1 ;
    }

    return clusters;
}

// Type clustering - 2 clusters with all fragments with only I and only D or 1 cluster with both (if we have at least one mixed fragment)
vector<Cluster> Insdeller::type_cluster(const Cluster &c) {
    vector<Cluster> tc;
    bool has_both = false;
    for (const Fragment &f : c) {
        has_both |= (f.t == 2) ;
    }

    if (has_both) {
        tc.push_back(c);
    } else {
        tc.push_back(Cluster(chrom)) ;
        tc.push_back(Cluster(chrom)) ;
        for (const Fragment &f : c) {
            Cluster &c = f.t == 0 ? tc.front() : tc.back();
            c.add_fragment(f);
        }
    }
    return tc ;
}

// string_view may help make this faster
unordered_map<string, int> get_unique_anchors(char* ref, int l, int k) {
    unordered_map<string, int> kmers ;
    unordered_map<string, int> positions ;
    char* _kmer = (char*) malloc(sizeof(char) * (k + 1)) ;
    _kmer[k] = '\0' ;
    for (int i = 0; i < l; i++) {
        memcpy(_kmer, ref + i, k) ;
        string kmer(_kmer) ;
        kmers[canonicalize(string(kmer))] += 1 ;
        positions[canonicalize(string(kmer))] = i ;
    }
    free(_kmer) ;
    auto it = kmers.begin() ;
    int n = 0 ;
    while (it != kmers.end()) {
        if (it->second != 1) {
           it = kmers.erase(it) ; 
           n += 1 ;
        } else {
            it++ ;
        }
    }
    return positions ;
}

Cluster Insdeller::merge_close_fragments(const Cluster& cluster, int distance) {
    Cluster merged_cluster(cluster.chrom) ;
    merged_cluster.s = cluster.s ;
    merged_cluster.e = cluster.e ;
    //
    unordered_map<string, vector<Fragment>> reads ;
    for (auto& f: cluster.fragments) {
        reads[f.name].push_back(f) ;
    }
    for (auto& r: reads) {
        auto& fragments = r.second ;
        for (auto i = fragments.begin(); i != fragments.end(); i++) {
            for (auto j = i++; j != fragments.end(); j++) {
                if (j->read_s - i->read_e < distance) {
                    i->ref_e = j->ref_e ;
                    i->read_e = j-> read_e ;
                    j = fragments.erase(j) ;
                    // TODO needs to update sequence
                }
            }
        }
        for (auto& f: fragments) {
            merged_cluster.add_fragment(f) ;
        }
    }
    return merged_cluster ;
}

bool Insdeller::should_filter_read(bam1_t* alignment, char* read_seq, string chrom, int* global_num_bases, int* global_num_mismatch) {
    auto cigar_offsets = decode_cigar(alignment) ;
    int l = 0 ;
    for (auto p: cigar_offsets) {
        l += p.first ;
    }
    //
    int n = 0 ;
    int m = 0 ;
    int ref_offset = alignment->core.pos ;
    int ins_offset = 0 ;
    int del_offset = 0 ;
    int match_offset = 0 ;
    int soft_clip_offset = 0 ;
    // Modify current bam1_t* struct
    int m_diff = 0 ;
    double num_match = 0 ;
    double num_mismatch = 0 ;
    while (true) {
        if (m == cigar_offsets.size()) {
            break ;
        }
        if (cigar_offsets[m].second == BAM_CMATCH || cigar_offsets[m].second == BAM_CEQUAL || cigar_offsets[m].second == BAM_CDIFF) {
            for (int j = 0; j < cigar_offsets[m].first; j++) {
                num_mismatch += 1 ? chromosome_seqs[chrom][ref_offset + j] != read_seq[match_offset + ins_offset + soft_clip_offset + j] : 0 ;
                n++ ;
            }
            ref_offset += cigar_offsets[m].first ;
            match_offset += cigar_offsets[m].first ;
            num_match += cigar_offsets[m].first ;
            m_diff = 0 ;
        } else if (cigar_offsets[m].second == BAM_CINS) {
            if (cigar_offsets[m].first <= 10) {
                // if a short INDEL then just don't add it to read
            } else {
                // for long INS, this is probably a SV so add it to the read
                for (int j = 0; j < cigar_offsets[m].first; j++) {
                    n++ ;
                }
            }
            ins_offset += cigar_offsets[m].first ;
        } else if (cigar_offsets[m].second == BAM_CDEL) {
            if (cigar_offsets[m].first <= 10) {
                // if a short DEL so let's just fix it
                for (int j = 0; j < cigar_offsets[m].first; j++) {
                    n++ ;
                }
                m_diff += cigar_offsets[m].first ;
            } else {
                // for long DEL, this is probably a SV so let it be what it was
            }
            del_offset += cigar_offsets[m].first ;
            ref_offset += cigar_offsets[m].first ;
        } else if (cigar_offsets[m].second == BAM_CSOFT_CLIP) {
            for (int j = 0; j < cigar_offsets[m].first; j++) {
                n++ ;
            }
            soft_clip_offset += cigar_offsets[m].first ;
        } else {
            // pass
        }
        m += 1 ;
    }
    *global_num_bases += num_match ;
    *global_num_mismatch += num_mismatch ;
    if (num_mismatch / num_match > 3 * expected_mismatch_rate) {
        return true ;
    }
    if (ins_offset + del_offset > 0.7 * strlen(read_seq)) {
        return true ;
    }
    return false ;
}

Cluster Insdeller::extend(Cluster &cluster, int extension) {
    std::sort(cluster.fragments.begin(), cluster.fragments.end()) ;
    auto c = cluster ;
    // find anchors
    int padding = 3000 ;
    int k = ANCHOR_SIZE ;
    char* ref = chromosome_seqs[c.chrom] + c.s - padding ;
    unordered_map<string, int> kmers ; //= 
    if (extension == EXTENSION_TYPE_KMER) { 
        cout << "Searching for anchors in " << c.chrom << "[" << c.s - padding << "-" << c.e + padding << "]" << endl ; 
        kmers = get_unique_anchors(ref, c.e - c.s + 2 * padding, k) ; 
        if (!kmers.size()) {
            lprint({"No anchors found for", c.get_id()}, 'E') ; 
            return cluster ;
        }
        cout << "Found" << kmers.size() << " anchors." << endl ;
    }

    int global_num_bases = 0 ;
    int global_num_mismatch = 0 ;

    unordered_map<string, std::vector<Fragment>> reads_in_cluster; // these are the reads we are interested in
    for (const Fragment& f : c) {
        reads_in_cluster[f.name].push_back(f) ;
    }
    //
    Cluster ext_c(chrom) ;
    string region = chrom + ":" + to_string(c.s) + "-" + to_string(c.e);
    hts_itr_t *itr = sam_itr_querys(read_bamindex[omp_get_thread_num()], read_bamhdr[omp_get_thread_num()], region.c_str()) ;
    bam1_t *aln = bam_init1() ;
    int buffer_len = 10000 ;
    uint8_t* seq = (uint8_t*) malloc(buffer_len) ;
    uint8_t* qual = (uint8_t*) malloc(buffer_len) ;
    char* _kmer = (char*) malloc(sizeof(char) * (k + 1)) ;
    _kmer[k] = '\0' ;
    int f = 0 ;
    int n = 0 ;
    while (sam_itr_next(read_bam[omp_get_thread_num()], itr, aln) > 0) {
        string qname(bam_get_qname(aln)) ;
        if (reads_in_cluster.find(qname) == reads_in_cluster.end()) {
            continue ;
        }
        // get read sequence
        uint32_t l = aln->core.l_qseq ; //length of the read
        if (l >= buffer_len) {
            free(seq) ;
            free(qual) ;
            buffer_len = l + 1 ;
            seq = (uint8_t*) malloc(buffer_len) ;
            qual = (uint8_t*) malloc(buffer_len) ;
        }
        uint8_t *s = bam_get_seq(aln) ; //quality string
        uint8_t *q = bam_get_qual(aln) ;
        for (int i = 0; i < l; i++) {
            seq[i] = seq_nt16_str[bam_seqi(s, i)] ; //gets nucleotide id and converts them into IUPAC id.
            qual[i] = 'I' ; //q[i] ;
        }
        seq[l] = '\0' ;
        qual[l] = '\0' ;
        vector<pair<int, int>> alpairs = get_aligned_pairs(aln) ;
        n++ ;
        if (extension == EXTENSION_TYPE_FULL && should_filter_read(aln, (char*) seq, cluster.chrom, &global_num_bases, &global_num_mismatch)) {
            f++ ;
            continue ;
        }
        // Extend each fragment on this read
        for (auto& f: reads_in_cluster[qname]) {
            // extend SFS until we reach anchors
            int new_ref_s ;
            int new_ref_e ;
            int new_read_s = -1 ;
            int new_read_e = -1 ;
            if (extension == EXTENSION_TYPE_KMER) { // kmer extension
                for (int i = f.read_s; i >= k - 1; i--) {
                    memcpy(_kmer, seq + i - (k - 1), k) ;
                    string canon = canonicalize(string(_kmer)) ;
                    if (kmers.find(canon) != kmers.end()) {
                        new_read_s = i - (k - 1) ;
                        break ;
                    }
                } 
                for (int i = f.read_e; i < l; i++) {
                    memcpy(_kmer, seq + i, k) ;
                    string canon = canonicalize(string(_kmer)) ;
                    if (kmers.find(canon) != kmers.end()) {
                        new_read_e = i + k ;
                        break ;
                    }
                }
                //cout << "Old: [" << f.read_s << " - " << f.read_e << "], [" << f.ref_s << ", " << f.ref_e << "]" << endl ;
                new_read_s = new_read_s == -1 ? f.read_s : new_read_s ;
                new_read_e = new_read_e == -1 ? f.read_e : new_read_e ;
                new_ref_s = f.ref_s - (f.read_s - new_read_s) ;
                new_ref_e = f.ref_e + (new_read_e - f.read_e) ;
            } else if (extension == EXTENSION_TYPE_FULL) { // full read
                new_read_s = 0 ;
                new_read_e = l - 1 ;
                new_ref_s = f.ref_s ;
                new_ref_e = f.ref_e ;
            } else if (extension == EXTENSION_TYPE_NONE) {
                new_read_s = f.read_s ;
                new_read_e = f.read_e ;
                new_ref_s = f.ref_s ;
                new_ref_e = f.ref_e ;
            } else {
                new_read_s = max(int(f.read_s) - extension, 0) ;
                new_read_e = min(int(l - 1), int(f.read_e) + extension) ;
                new_ref_s = f.ref_s - extension ;
                new_ref_e = f.ref_e + extension ;
            }
            string subseq((char*) seq, new_read_s, new_read_e - new_read_s + 1);
            string subqual((char*) qual, new_read_s, new_read_e - new_read_s + 1);
            Fragment new_f(qname, new_ref_s, new_ref_e, new_read_s, new_read_e, 0, subseq, subqual) ;
            ext_c.add_fragment(new_f) ;
        }
    }
    cout << "Filtered " << f << " fragments out of " << n << ". " << ext_c.size() << " remaining." << endl ;
    free(seq) ;
    free(qual) ;
    free(_kmer) ;
    return ext_c ;
}

int Insdeller::cluster_anchors(const Cluster &cluster) {
    unordered_map<string, vector<int>> left ;
    unordered_map<string, vector<int>> right ;
    for (int i = 0; i < cluster.fragments.size(); i++) {
        auto& f = cluster.fragments[i] ;
        string r ;
        string l ;
        if (f.right_extended) {
            r = f.seq.substr(0, ANCHOR_SIZE) ;
            right[r].push_back(i) ;
        }
        if (f.left_extended) {
            l = f.seq.substr(f.seq.length() - ANCHOR_SIZE, ANCHOR_SIZE) ;
            left[l].push_back(i) ;
        }
    }
    if (left.size() == 1 && right.size() == 1) {
        // All have been extended to the same anchor, easy
        return 0 ;
    } else if (left.size() == 1 || right.size() == 1) {
        // still easy
        return 1 ;
    } else {
        return max(left.size(), right.size()) ;
    }
}

// puts overlapping fragments of the same size in the same cluster
vector<Cluster> Insdeller::cluster_breakpoints(const Cluster& cluster, float ratio) {
    //
    ratio = 0.7 ;
    vector<Cluster> clusters ;
    // merge overlapping fragments of similar sizes
    for (const auto& f: cluster) {
        bool new_cluster = true ;
        vector<int> to_merge ;
        for (int i = 0; i < clusters.size(); i++) {
            int ref_s = clusters[i].s ;
            int ref_e = clusters[i].e ;
            int l_2 = ref_e - ref_s ;
            int l_1 = f.ref_e - f.ref_s ;
            if (min(int(f.ref_e), ref_e) - max(int(f.ref_s), ref_s) >= 0) {
                if (min(l_1, l_2) / max(l_1, l_2) >= ratio) {
                    to_merge.push_back(i) ;
                    new_cluster = false ;
                }
            }
        }
        if (new_cluster) {
            Cluster c(cluster.chrom) ;
            c.add_fragment(f) ;
            clusters.push_back(c) ;
        } else {
            int i = to_merge[0] ;
            int ref_s = clusters[i].s ;
            int ref_e = clusters[i].e ;
            for (int j = 1; j < to_merge.size(); j++) {
                int index = to_merge[j] ;
                for (const auto& ff: clusters[index].fragments) {
                    clusters[i].add_fragment(ff) ;
                }
                clusters[index].fragments.clear() ;
            }
            clusters[i].add_fragment(f) ;
            int offset = 0 ;
            for (int j = 1; j < to_merge.size(); j++) {
                clusters.erase(clusters.begin() + to_merge[j] - offset) ;
                offset += 1 ;
            }
        }
    }
    // merge clusters with each other
    int last_b = -1 ;
    int to_remove[clusters.size()] ;
    for (int i = 0; i < clusters.size(); i++) {
        to_remove[i] = 0 ;
    }
    while (last_b != clusters.size()) {
        std::sort(clusters.begin(), clusters.end()) ;
        last_b = clusters.size() ;
        for (int _ = 0; _ != clusters.size(); _++) {
            if (to_remove[_] == 1) {
                continue ;
            }
            auto& c = clusters[_] ;
            vector<int> to_merge ;
            for (int i = 0; i < clusters.size(); i++) {
                if (clusters[i] == c || to_remove[i] == 1) {
                    continue ;
                }
                int s = clusters[i].s ;
                int e = clusters[i].e ;
                int l_1 = c.e - c.s ;
                int l_2 = e - s ;
                int o = min(int(c.e), e) - max(int(c.s), s) ;
                float r = float(min(l_1, l_2)) / float(max(l_1, l_2)) ;
                if (o >= 0) {
                    //if ((l_1 <= 10 && l_2 <= 10) || r >= ratio || o == min(l_1, l_2) || o >= ratio * min(l_1, l_2)) {
                        to_merge.push_back(i) ;
                        //cout << "Comparing " << c.get_id() << " to " << clusters[i].get_id() << " o: " << o << ", r: " << r << endl ;
                    //}
                }
            }
            if (to_merge.size() != 0) {
                for (int j = 0; j < to_merge.size(); j++) {
                    int index = to_merge[j] ;
                    for (const auto& ff: clusters[index].fragments) {
                        c.add_fragment(ff) ;
                    }
                    to_remove[index] = 1 ;
                    clusters[index].fragments.clear() ;
                }
            }
        }
    }
    int i = 0 ;
    auto it = clusters.begin() ;
    while (it != clusters.end()) {
        if (it->fragments.size() < 2 || to_remove[i] == 1) {
            it = clusters.erase(it) ;
        } else {
            it++ ;
        }
        i++ ;
    }
    cout << "Clustering " << cluster.fragments.size() << " fragments in " << cluster.get_id() << " into " << clusters.size() << " clusters." << endl ;
    //for (auto& c: clusters) {
    //    cout << c.chrom << ":" << c.s << "-" << c.e << endl ;
    //}
    return clusters ;
}

vector<Cluster> Insdeller::scluster(const Cluster &cluster) {
    vector<Cluster> clusters ;

    // greedy choice: longer fragment is representative #1
    Fragment repr1;
    for (const Fragment &f : cluster) {
        if (f.size() > repr1.size()) {
            repr1 = f;
        }
    }

    // the other representative will be the fragment with the lowest similarity with repr1
    Fragment repr2 ;
    double min_d = 100.0 ;
    for (const Fragment &f : cluster) {
        double d = rapidfuzz::fuzz::ratio(repr1.seq, f.seq);
        if (d < 60 && d < min_d) {
            min_d = d;
            repr2 = f;
        }
    }

    if (repr2.size() == 0) {
        clusters.push_back(cluster);
    }
    else {
        clusters.push_back(Cluster(cluster.chrom));
        clusters.push_back(Cluster(cluster.chrom));
        double d1, d2;
        for (const Fragment &f : cluster) {
            d1 = rapidfuzz::fuzz::ratio(repr1.seq, f.seq);
            d2 = rapidfuzz::fuzz::ratio(repr2.seq, f.seq);
            if (d1 <= d2) {
                clusters.front().add_fragment(f);
            } else {
                clusters.back().add_fragment(f);
            }
        }
    }
    return clusters;
}


CIGAR Insdeller::align_ksw2(const char *ref, const char *query, int sc_mch, int sc_mis, int gapo, int gape) {
    int i, a = sc_mch, b = sc_mis < 0 ? sc_mis : -sc_mis; // a>0 and b<0
    int8_t mat[25] = {a, b, b, b, 0, b, a, b, b, 0, b, b, a, b, 0, b, b, b, a, 0, 0, 0, 0, 0, 0};
    int tl = strlen(ref), ql = strlen(query);
    uint8_t *ts, *qs;
    ksw_extz_t ez;

    memset(&ez, 0, sizeof(ksw_extz_t)) ;
    ts = (uint8_t*) malloc(tl) ;
    qs = (uint8_t*) malloc(ql) ;
    for (i = 0; i < tl; ++i) {
        ts[i] = _nt4_table[(uint8_t) ref[i]] ; // encode to 0/1/2/3
    }
    for (i = 0; i < ql; ++i) {
        qs[i] = _nt4_table[(uint8_t) query[i]] ;
    }
    // ksw_extz(0, ql, qs, tl, ts, 5, mat, gapo, gape, -1, -1, 0, &ez);

    ksw_extd(0, ql, qs, tl, ts, 5, mat, 50, 1, 50, 0, -1, -1, 0, &ez) ;

    vector<pair<uint, char>> cigar(ez.n_cigar);
    for (i = 0; i < ez.n_cigar; ++i) // print CIGAR
        cigar[i] = make_pair(ez.cigar[i] >> 4, "MID"[ez.cigar[i] & 0xf]);
    int score = ez.score ;
    free(ez.cigar) ;
    free(ts) ;
    free(qs) ;
    return CIGAR(cigar, score, 0) ;
}

CIGAR Insdeller::align_edlib(const char *ref, const char *query, int sc_mch, int sc_mis, int gapo, int gape) {
    EdlibAlignResult result = edlibAlign(query, strlen(query), ref, strlen(ref),
        edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
    if (result.status == EDLIB_STATUS_OK) {
        char* cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
        CIGAR c(cigar, result.editDistance, result.startLocations[0]) ;
        free(cigar) ;
        edlibFreeAlignResult(result) ;
        return c ;
    } else {
        edlibFreeAlignResult(result) ;
        return CIGAR({}, -1, 0) ;
    }
}

Cluster Insdeller::compress_cluster(const Cluster& c, int size) {
    if (c.size() <= size) {
        return c ;
    }
    cout << "Compressing cluster with " << c.size() << " fragments.." << endl ;
    Cluster cluster(c.chrom) ;
    unordered_map<int, bool> indices ;
    while (true) {
        int r = rand() % c.size() ;
        if (indices.find(r) == indices.end()) {
            indices[r] = true ;
            cluster.add_fragment(c.fragments[r]) ;
        } else {
            continue ;
        }
        if (cluster.fragments.size() == size) {
            break ;
        }
    }
    cout << cluster.fragments.size() << " remaining." << endl ;
    return cluster ;
}

vector<SV> Insdeller::call_poa_svs(Cluster &cluster, const string &ref, ofstream &o) {
    int padding = 0 ;
    vector<string> _consensus = cluster.poa() ;
    vector<SV> svs ;
    cout << _consensus.size() << " POA assemblies." << endl ;
    for (auto consensus: _consensus) {
        cout << consensus << endl ;
    }
    for (auto consensus: _consensus) {
        vector<SV> _svs = Haplotyper().map_to_chm13(cluster, consensus) ;
        if (_svs.size() != 0) {
            cout << "[E]" << " found SVs in POA mapping to CHM13. Consensus is invalid." << endl ;
            for (const auto& sv: _svs) {
                cout << sv << endl ;
            }
            int a ; 
            cin >> a ;
        }
        padding = max(0, int(consensus.length()) - int(cluster.e - cluster.s)) ;
        string reference = string(ref, cluster.s - padding, cluster.e - cluster.s + 2 * padding) ; 
        cout << "Reference length:" << reference.length() << ", Consensus length: " << consensus.length() << endl ; 
        CIGAR cigar = align_ksw2(reference.c_str(), consensus.c_str(), 1, -3, 100, 0); // TODO: adjust scores/penalties. I want less gaps
        if (cigar.score == -1) {
            cout << "Local alignment fail, trying global alignemnt.." << endl ;
        }
        cigar.fixclips() ;
        cigar.print() ;
    
        o << chrom << ":" << cluster.s + 1 << "-" << cluster.e + 1 << ":" << cluster.size() << "\t"
            << 0 << "\t"
            << cluster.chrom << "\t"
            << cluster.s + 1 << "\t"
            << 60 << "\t"
            << cigar.to_str() << "\t"
            << "*"
            << "\t"
            << 0 << "\t"
            << 0 << "\t"
            << consensus << "\t"
            << "*"
            << endl;

        uint qs = 0 ;
        uint rs = cluster.s - padding + cigar.start ;
        bool saw_match = false ;
        for (uint i = 0; i < cigar.size(); i++) {
            SV sv ;
            if (cigar[i].second == 'S') {
                qs += cigar[i].first ;
                saw_match = true ;
            } else if (cigar[i].second == 'I') {
                string ref_allele = ref.substr(rs - 1, 1) ;
                string alt_allele = ref_allele + consensus.substr(qs, cigar[i].first) ;
                if (saw_match) {
                    sv = SV("INS", cluster.chrom, rs - 1, ref_allele, alt_allele, cluster.size(), cluster.full_cov, cigar.ngaps, cigar.score) ;
                }
                qs += cigar[i].first ;
            } else if (cigar[i].second == 'D') {
                string ref_allele = ref.substr(rs - 1, cigar[i].first + 1) ;
                string alt_allele = ref_allele.substr(0, 1) ;
                if (saw_match) {
                    sv = SV("DEL", cluster.chrom, rs - 1, ref_allele, alt_allele, cluster.size(), cluster.full_cov, cigar.ngaps, cigar.score) ;
                }
                rs += cigar[i].first ;
            } else if (cigar[i].second == 'M') {
                saw_match = true ;
                rs += cigar[i].first ;
                qs += cigar[i].first ;
            } else {
                cerr << "Unknown CIGAR op " << cigar[i].second << endl;
                exit(1);
            }
            if (abs(sv.l) > 25 && (sv.s >= cluster.s - 100 && sv.e <= cluster.e + 100)) {
                cout << sv << endl ;
                svs.push_back(sv) ;
            }
        }
        cout << endl ;
    }
    sort(svs.begin(), svs.end()) ;
    int a ; 
    return svs ;
}

vector<SV> Insdeller::call_svs(const Cluster& cluster, const string& ref) {
    cout << cluster.get_id() << endl ;
    int coverage = 0 ;
    //bam1_t *aln = bam_init1() ;
    //string region = cluster.chrom + ":" + to_string(cluster.s) + "-" + to_string(cluster.e + 1) ;
    //hts_itr_t *itr = sam_itr_querys(read_bamindex[omp_get_thread_num()], read_bamhdr[omp_get_thread_num()], region.c_str()) ;
    //while (sam_itr_next(read_bam[omp_get_thread_num()], itr, aln) > 0) {
    //    if (aln->core.flag & BAM_FUNMAP || aln->core.flag & BAM_FSUPPLEMENTARY || aln->core.flag & BAM_FSECONDARY) {
    //        continue ;
    //    }
    //    coverage++ ;
    //}
    int distance = 10 ;
    vector<SV> _svs ;
    for (const Fragment &rep: cluster) {
        int last_r = -1 ;
        int last_q = -1 ;
        bool found_del = true ;
        for (int i = 0; i < rep.aligned_pairs.size(); i++) {
            int q = rep.aligned_pairs[i].first ;
            int r = rep.aligned_pairs[i].second ;
            if (q == -1) {
                found_del = true ;
                continue ;
            } else {
                if (found_del) {
                    found_del = false ;
                    if (r != -1 && last_r != -1 && r != last_r + 1) {
                        int svlen = r - last_r - 1 ;
                        string ref_allele = ref.substr(last_r, svlen) ;
                        string alt_allele = ref_allele.substr(0, 1) ;
                        SV sv = SV("DEL", cluster.chrom, last_r, ref_allele, alt_allele, cluster.size(), coverage, -1, -1) ;
                        if (svlen > 10) {
                            _svs.push_back(sv) ;
                        }
                    }
                }
            }
            last_q = q ;
            last_r = r ; 
        }
        // find insertions
        last_q = -1 ;
        last_r = -1 ;
        bool found_ins = true ;
        for (int i = 0; i < rep.aligned_pairs.size(); i++) {
            int q = rep.aligned_pairs[i].first ;
            int r = rep.aligned_pairs[i].second ;
            if (r == -1) {
                 found_ins = true ;
                 continue ;
            } else {
                if (found_ins) {
                    found_ins = false ;
                    if (q != -1 && last_q != -1 && q != last_q + 1) {
                        int svlen = q - last_q - 1 ;
                        string ref_allele = ref.substr(last_r, 1) ;
                        string alt_allele = ref_allele + rep.seq.substr(last_q + 1, svlen) ;
                        SV sv = SV("INS", cluster.chrom, last_r, ref_allele, alt_allele, cluster.size(), coverage, -1, -1) ;
                        if (svlen > 10) {
                            _svs.push_back(sv) ;
                        }
                    }
                }
            }
            last_q = q ;
            last_r = r ; 
        }
    }
    // cluster the SV, ignore SVs predicted by less than 2 fragments, merge close SV with similar size to the more common SV
    sort(_svs.begin(), _svs.end()) ;
    vector<SVCluster> sv_clusters ;
    for (int _ = 0; _ < _svs.size(); _++) {
        auto sv = _svs[_] ;
        //cout << sv << endl ;
        vector<int> to_merge ;
        bool new_cluster = true ;
        for (int j = 0; j < sv_clusters.size(); j++) {
            auto& sv_c = sv_clusters[j] ;
            int l_1 = abs(sv.l) ;
            int l_2 = sv_c.svlen ;
            if (sv_c.type == sv.type) {
                if (abs(sv_c.s - int(sv.s)) <= distance) {
                    if (min(l_1, l_2) / max(l_1, l_2) >= 0.9) {
                        new_cluster = false ;
                        to_merge.push_back(j) ;
                    }
                }
            }
        }
        if (new_cluster) {
            sv_clusters.push_back(SVCluster(sv)) ;
        } else {
            int i = to_merge[0] ;
            for (int j = 1; j < to_merge.size(); j++) {
                int index = to_merge[j] ;
                for (auto& sv: sv_clusters[index].svs) {
                    sv_clusters[i].add_sv(sv.first, sv.second) ;
                }
            }
            sv_clusters[i].add_sv(sv) ;
            int offset = 0 ;
            for (int j = 1; j < to_merge.size(); j++) {
                sv_clusters.erase(sv_clusters.begin() + to_merge[j] - offset) ;
                offset += 1 ;
            }
        }
    }
    int last_size = -1 ;
    int to_remove[sv_clusters.size()] ;
    memset(&to_remove, 0, sv_clusters.size() * sizeof(int)) ;
    while (last_size != sv_clusters.size()) {
        last_size = sv_clusters.size() ;
        for (int j = 0; j != sv_clusters.size(); j++) {
            if (to_remove[j] == 1) {
                continue ;
            }
            auto& sv_c = sv_clusters[j] ;
            vector<int> to_merge ;
            for (int i = 0; i < sv_clusters.size(); i++) {
                if (i == j || to_remove[i] == 1 || sv_clusters[i].type != sv_c.type) {
                    continue ;
                }
                int s = sv_clusters[i].s ;
                int l_1 = sv_c.svlen ;
                int l_2 = sv_clusters[i].svlen ;
                float r = float(min(l_1, l_2)) / float(max(l_1, l_2)) ;
                if (abs(sv_clusters[i].s - sv_c.s) <= distance) {
                    to_merge.push_back(i) ;
                }
            }
            if (to_merge.size() != 0) {
                for (int j = 0; j < to_merge.size(); j++) {
                    int index = to_merge[j] ;
                    for (auto& sv: sv_clusters[index].svs) {
                        sv_c.add_sv(sv.first, sv.second) ;
                    }
                    to_remove[index] = 1 ;
                }
            }
        }
    }
    vector<SV> svs ;
    for (int i = 0; i < sv_clusters.size(); i++) {
        if (to_remove[i] == 1) {
            continue ;
        }
        auto sv_c = sv_clusters[i] ;
        int m = 0 ;
        SV sv ;
        for (auto& kv: sv_c.svs) {
            if (kv.second > m) {
                sv = kv.first ;
                m = kv.second ;
            }
        }
        if (m >= 3) {
            //cout << sv << endl ;
            sv.w = m ;
            svs.push_back(sv) ;
        }
    }
    cout << "Predicted " << svs.size() << " SVs." << endl ;
    return svs ;
}

vector<SV> intersect_svs(vector<SV> a, vector<SV> b) {
    vector<SV> svs ;
    for (auto& s: a) {
        for (int j = 0; j < b.size(); j++) {
            SV t = b[j] ;
            cout << "Matching " << s.idx << " with " << t.idx << endl ;
            if (abs(s.s - t.s) < 10) {
                int l_1 = abs(s.l) ;
                int l_2 = abs(t.l) ;
                float r = float(min(l_1, l_2)) / float(max(l_1, l_2)) ;
                if (r >= 0.95 || abs(l_1 - l_2) <= 5) {
                    svs.push_back(s) ;
                    cout << "Matched " << s.idx << " with " << t.idx << endl ;
                    b.erase(b.begin() + j) ;
                    break ;
                }
            }
        }
    }
    cout << "Intersecting " << a.size() << " SVs with " << b.size() << " SVs; " << svs.size() << " remain." << endl ;
    return svs ;
}

vector<SV> Insdeller::filter_chain_svs(vector<SV> svs) {
    sort(svs.begin(), svs.end()) ;
    vector<SV> deletions ;
    vector<SV> insertions ;
    for (int i = 0; i < svs.size(); i++) {
        if (svs[i].type == "DEL") {
            deletions.push_back(svs[i]) ;
        } else {
            insertions.push_back(svs[i]) ;
        }
    }
    int i = 0 ;
    cout << svs.size() << " SV before chain filtering.." << endl ;
    vector<SV> final_svs ;
    for (auto& _svs: {deletions, insertions}) {
        while (i < _svs.size()) {
            cout << "Looking for chain starting at " << _svs[i].idx << endl ;
            // this will find a chain of SVs of similar size
            // e.g see chr21:8,459,874-8,462,286 on CHM13 mapped to hg38
            vector<SV> sv_cluster ;
            sv_cluster.push_back(_svs[i]) ;
            int j = i + 1 ;
            while (j < _svs.size()) {
                cout << "Extending to " << _svs[j].idx << endl ;
                int l_1 = abs(sv_cluster.back().l) ;
                int l_2 = abs(_svs[j].l) ;
                bool found = false ;
                int d = abs(_svs[j].s - sv_cluster.back().s) ;
                float r = float(min(l_1, l_2)) / float(max(l_1, l_2)) ;
                cout << d << " " << r << endl ;
                if (d < 500) {
                    if (r > 0.95) {
                        sv_cluster.push_back(_svs[j]) ;
                        found = true ;
                    }
                } else {
                    if (r > 0.99) {
                        // We don't expect to see two SVs of exact same size after each other
                        // This seems to improve things
                        sv_cluster.push_back(_svs[j]) ;
                        found = true ;
                    }
                }
                if (!found) {
                    break ;
                }
                j++ ;
            }
            i = j ;
            cout << "Found SV chain of svlen " << sv_cluster[0].l << " with " << sv_cluster.size() << " SVs." << endl ;
            int m = 0 ;
            SV sv ;
            for (SV s: sv_cluster) {
                if (s.w > m) {
                    m = s.w ;
                    sv = s ;
                }
            }
            final_svs.push_back(sv) ;
            cout << sv << endl ;
        }
    }
    sort(final_svs.begin(), final_svs.end()) ;
    // heuristic
    // Maybe do another round when we filter consecutive SVs of similar size
    cout << final_svs.size() << " SVs remain." << endl ;
    return final_svs ;
}

vector<SV> Insdeller::remove_duplicate_svs(const vector<SV> &svs) {
    vector<SV> usvs ;
    set<string> uidxs ;
    for (const SV &sv: svs) {
        if (uidxs.find(sv.idx) == uidxs.end()) {
            uidxs.insert(sv.idx) ;
            usvs.push_back(sv) ;
        }
    }
    return usvs ;
}

vector<SV> Insdeller::call_batch(vector<Cluster>& position_clusters, const string& chrom_seq, ofstream& osam) {
    vector<SV> svs ;
    for (auto& pc: position_clusters) {
        if (pc.size() < 2) {
            continue ;
        }
        vector<Cluster> type_clusters = type_cluster(pc);
        for (Cluster &tc : type_clusters) {
            //if (tc.size() < 2) {
            //    continue ;
            //}
            //if (tc.s > 5300000) { //|| tc.e > 5280000) {
            //    continue ;
            //}
            //if (tc.get_id() != "chr21:7933252-7946087") {
            //    continue ;
            //}
            auto breakpoints = cluster_breakpoints(tc, 0.7) ;
            if (breakpoints.size() == 0) {
                tc.clear() ;
                break ;
            }
            if (breakpoints.size() == 1 && tc.get_type() != 2) {
                cout << omp_get_thread_num() << " " << "Using normal genotyper for " << tc.get_id() << endl ;
                auto cluster_svs = call_svs(breakpoints[0], chrom_seq) ;
                svs.insert(svs.begin(), cluster_svs.begin(), cluster_svs.end()) ;
            } else if (breakpoints.size() < 25) {
                cout << omp_get_thread_num() << "Delegating to assembler for " << tc.get_id() << endl ;
                // 1. SFS-only POA
                //for (auto breakpoint: breakpoints) {
                //    Cluster extended_breakpoint = extend(breakpoint, EXTENSION_TYPE_KMER) ;
                //    cout << extended_breakpoint.get_id() << endl ;
                //    auto _poa_svs = call_poa_svs(breakpoint, chrom_seq, osam) ;
                //    svs.insert(svs.begin(), _poa_svs.begin(), _poa_svs.end()) ;
                //}
                // 2. Full-read POA
                //auto full_extended_tc = extend(tc, EXTENSION_TYPE_FULL) ;
                //auto compressed_tc = compress_cluster(full_extended_tc, 100) ;
                //if (compressed_tc.size() < 2) {
                //    continue ;
                //}
                //auto _poa_svs = call_poa_svs(compressed_tc, chrom_seq, osam) ;
                //svs.insert(svs.begin(), _poa_svs.begin(), _poa_svs.end()) ;
                // 3. Full-read miniasm
                auto full_extended_tc = extend(tc, EXTENSION_TYPE_FULL) ;
                auto _miniasm_svs = Haplotyper().haplotype(full_extended_tc) ;
                svs.insert(svs.begin(), _miniasm_svs.begin(), _miniasm_svs.end()) ;
            }
            breakpoints.clear() ;
            tc.clear() ;
        }
        pc.clear() ;
    }
    return svs ;
}

void Insdeller::call(const string &chrom_seq, ofstream &osam) {
    srand(9216429) ;
    vector<Cluster> _position_clusters = position_cluster() ;
    cout << "Sorted reads into " << _position_clusters.size() << " P-clusters." << endl ;
    vector<vector<Cluster>> position_clusters(config->threads) ;
    int t = 0 ;
    for (auto& c: _position_clusters) {
        position_clusters[t % config->threads].push_back(c) ;
        t++ ;
    }
    _position_clusters.clear() ;

    vector<vector<SV>> svs(config->threads) ;
    int m = 0 ;
    #pragma omp parallel for num_threads(config->threads)
    for (int i = 0; i < config->threads; i++) {
        svs[i] = call_batch(position_clusters[i], chrom_seq, osam) ;
    }
    cout << "Done.." << endl ;
    vector<SV> merged_svs ; 
    for (int i = 0; i < config->threads; i++) {
        for (const SV &sv : svs[i]) {
            merged_svs.push_back(sv) ;
        }
    }
    osvs = filter_chain_svs(remove_duplicate_svs(merged_svs)) ;
    cout << "Discovered " << osvs.size() << " SVs.. " << endl ;
    for (int i = 0; i < config->threads; i++) {
        sam_close(read_bam[i]) ;
    }
}
