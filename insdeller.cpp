#include "insdeller.hpp"

#define ANCHOR_SIZE 7

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
            uint op = bam_cigar_op(*(cigar + i));
            has_i |= op == BAM_CINS;
            has_d |= op == BAM_CDEL;
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


Cluster Insdeller::extend(Cluster &cluster) {
    // sort cluster fragments based on position
    std::sort(cluster.fragments.begin(), cluster.fragments.end()) ;
    // TODO: merge close signatures together?
    // auto c = merge_close_fragments(cluster, 50) ;
    auto c = cluster ;
    // find anchors
    int padding = 3000 ;
    int k = ANCHOR_SIZE ;
    char* ref = chromosome_seqs[c.chrom] + c.s - padding ;
    //cout << "Searching for anchors in " << c.chrom << "[" << c.s - padding << "-" << c.e + padding << "]" << endl ; 
    auto kmers = get_unique_anchors(ref, c.e - c.s + 2 * padding, k) ;
    //cout << "Found " << kmers.size() << " anchors.." << endl ;
    if (!kmers.size()) {
        lprint({"No anchors found for", c.get_id()}, 'E') ; 
    }

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
    while (sam_itr_next(read_bam[omp_get_thread_num()], itr, aln) > 0) {
        bool left_extend = false ;
        bool right_extend = false ;
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
        // Extend each fragment on this read
        for (auto& f: reads_in_cluster[qname]) {
            // extend SFS until we reach anchors
            int new_read_s = -1 ;
            int new_read_e = -1 ;
            if (!config->assemble) {
                for (int i = f.read_s - 1; i >= k - 1; i--) {
                    memcpy(_kmer, seq + i - (k - 1), k) ;
                    string canon = canonicalize(string(_kmer)) ;
                    if (kmers.find(canon) != kmers.end()) {
                        new_read_s = i - (k - 1) ;
                        left_extend = true ;
                        break ;
                    }
                } 
                for (int i = f.read_e; i < l; i++) {
                    memcpy(_kmer, seq + i, k) ;
                    string canon = canonicalize(string(_kmer)) ;
                    if (kmers.find(canon) != kmers.end()) {
                        new_read_e = i + k ;
                        right_extend = true ;
                        break ;
                    }
                } 
                if (new_read_s == -1 || new_read_e == -1) {
                    //lprint({"Anchor extension unsuccessful for SFS", qname, to_string(f.read_s), "-", to_string(f.read_e)}, 'E') ;
                }
            }
            //cout << "Old: [" << f.read_s << " - " << f.read_e << "], [" << f.ref_s << ", " << f.ref_e << "]" << endl ;
            new_read_s = new_read_s == -1 ? f.read_s : new_read_s ;
            new_read_e = new_read_e == -1 ? f.read_e : new_read_e ;
            // TODO: for minasm only
            if (config->assemble) {
                new_read_s = 0 ;
                new_read_e = l - 1 ;
            }
            // update alignment
            int last_r = aln->core.pos ;
            int new_ref_s = -1 ;
            int new_ref_e = -1 ;
            // This will work if we haven't extended the kmer
            // If we have extended the kmer, then it may fail
            for (uint i = 0; i < alpairs.size(); i++) {
                int q = alpairs[i].first ;
                int r = alpairs[i].second ;
                if (r != -1) {
                    last_r = r ;
                }
                if (q == new_read_s) {
                    new_ref_s = r != -1 ? r : last_r ;
                } else if (q >= new_read_e) {
                    if (new_ref_e == -1) {
                        new_ref_e = r ;
                        if (r != -1) {
                            break ;
                        }
                    }
                }
            }
            if (new_read_e == l - 1) {
                new_ref_e = bam_endpos(aln) ;
            }
            if (new_ref_s == -1) {
                new_ref_s = f.ref_s ;
            }
            if (new_ref_e == -1) {
                new_ref_e = f.ref_e ;
            }
            if (new_ref_s == -1 && new_ref_e == -1) {
                cerr << "[W] no reference base on " << chrom << ":" << c.s << "-" << c.e << " (" << qname << ")" << endl;
                continue ;
            }
            //cout << "New: [" << new_read_s << " - " << new_read_e <<"], [" << new_ref_s << ", " << new_ref_e << "]" << endl ;
            //TODO why the heck is this happening?
            if (new_ref_s > f.ref_s || new_ref_e < f.ref_e) {
                //lprint({"Invalid SFS mapping for", f.name}, 'E') ;
                continue ;
            }
            //assert(new_ref_s <= f.ref_s) ;
            //assert(new_ref_e >= f.ref_e) ;
            string subseq((char*) seq, new_read_s, new_read_e - new_read_s + 1);
            string subqual((char*) qual, new_read_s, new_read_e - new_read_s + 1);
            Fragment new_f(qname, new_ref_s, new_ref_e, new_read_s, new_read_e, 0, subseq, subqual) ;
            new_f.o_ref_s = f.ref_s ;
            new_f.o_ref_e = f.ref_e ;
            new_f.extend(left_extend, right_extend) ;
            ext_c.add_fragment(new_f) ;
        }
    }
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
vector<Cluster> Insdeller::cluster_breakpoints(Cluster& cluster, float ratio) {
    //
    ratio = 0.7 ;
    vector<Cluster> clusters ;
    // merge overlapping fragments of similar sizes
    for (auto& f: cluster) {
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
                for (auto& ff: clusters[index].fragments) {
                    clusters[i].add_fragment(ff) ;
                }
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
    memset(&to_remove, 0, clusters.size() * sizeof(int)) ;
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
                    if ((l_1 <= 10 && l_2 <= 10) || r >= ratio) {
                        to_merge.push_back(i) ;
                        //cout << "Comparing " << c.get_id() << " to " << clusters[i].get_id() << " o: " << o << ", r: " << r << endl ;
                    }
                }
            }
            if (to_merge.size() != 0) {
                for (int j = 0; j < to_merge.size(); j++) {
                    int index = to_merge[j] ;
                    for (auto& ff: clusters[index].fragments) {
                        c.add_fragment(ff) ;
                    }
                    to_remove[index] = 1 ;
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
    Fragment repr2;
    double min_d = 100.0;
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

CIGAR Insdeller::align(const char *tseq, const char *qseq, int sc_mch, int sc_mis, int gapo, int gape) {
    int i, a = sc_mch, b = sc_mis < 0 ? sc_mis : -sc_mis; // a>0 and b<0
    int8_t mat[25] = {a, b, b, b, 0, b, a, b, b, 0, b, b, a, b, 0, b, b, b, a, 0, 0, 0, 0, 0, 0};
    int tl = strlen(tseq), ql = strlen(qseq);
    uint8_t *ts, *qs;
    ksw_extz_t ez;

    memset(&ez, 0, sizeof(ksw_extz_t)) ;
    ts = (uint8_t*) malloc(tl) ;
    qs = (uint8_t*) malloc(ql) ;
    for (i = 0; i < tl; ++i) {
        ts[i] = _nt4_table[(uint8_t) tseq[i]] ; // encode to 0/1/2/3
    }
    for (i = 0; i < ql; ++i) {
        qs[i] = _nt4_table[(uint8_t) qseq[i]] ;
    }
    // ksw_extz(0, ql, qs, tl, ts, 5, mat, gapo, gape, -1, -1, 0, &ez);

    ksw_extd(0, ql, qs, tl, ts, 5, mat, 50, 1, 50, 0, -1, -1, 0, &ez) ;

    vector<pair<uint, char>> cigar(ez.n_cigar) ;
    for (i = 0; i < ez.n_cigar; ++i) {// print CIGAR
        cigar[i] = make_pair(ez.cigar[i] >> 4, "MID"[ez.cigar[i] & 0xf]);
    }
    int score = ez.score ;
    free(ez.cigar) ;
    free(ts) ;
    free(qs) ;
    return CIGAR(cigar, score) ;
}

vector<SV> Insdeller::call_poa_svs(const Cluster &c, const string &ref, ofstream &o) {
    cout << "POA: " << c.get_id() << endl ;
    vector<SV> svs;
    string reference(ref, c.s, c.e - c.s + 1);
    string consensus = c.poa();
    cout << "Consensus length: " << consensus.length() << endl ; 
    CIGAR cigar = align(reference.c_str(), consensus.c_str(), 1, -3, 100, 0); // TODO: adjust scores/penalties. I want less gaps
    cigar.fixclips();

    // SAM output
    o << chrom << ":" << c.s + 1 << "-" << c.e + 1 << ":" << c.size() << "\t"
      << 0 << "\t"
      << c.chrom << "\t"
      << c.s + 1 << "\t"
      << 60 << "\t"
      << cigar.to_str() << "\t"
      << "*"
      << "\t"
      << 0 << "\t"
      << 0 << "\t"
      << consensus << "\t"
      << "*"
      << endl;
    // << "NM:i:" << 0 << endl;
    // << "AS:i:" << cigar.score << endl;

    // extracting insertions/deletions breakpoints
    uint qs = 0 ;       // position on consensus
    uint lrs = 0 ;      // position on local reference
    uint rs = c.s ; // position on reference

    for (uint i = 0; i < cigar.size(); i++) {
        SV sv ;
        if (cigar[i].second == 'S') {
            qs += cigar[i].first ;
        } else if (cigar[i].second == 'I') {
            string ref_allele = ref.substr(rs - 1, 1) ;
            string alt_allele = ref_allele + consensus.substr(qs, cigar[i].first) ;
            sv = SV("INS", c.chrom, rs - 1, ref_allele, alt_allele, c.size(), c.full_cov, cigar.ngaps, cigar.score) ;
            qs += cigar[i].first ;
        } else if (cigar[i].second == 'D') {
            string ref_allele = ref.substr(rs - 1, cigar[i].first + 1) ;
            string alt_allele = ref_allele.substr(0, 1) ;
            sv = SV("DEL", c.chrom, rs - 1, ref_allele, alt_allele, c.size(), c.full_cov, cigar.ngaps, cigar.score) ;
            rs += cigar[i].first ;
            lrs += cigar[i].first ;
        } else if (cigar[i].second == 'M') {
            rs += cigar[i].first ;
            qs += cigar[i].first ;
            lrs += cigar[i].first ;
        } else {
            cerr << "Unknown CIGAR op " << cigar[i].second << endl;
            exit(1);
        }
        if (abs(sv.l) > 0) {
            cout << sv << endl ;
            svs.push_back(sv) ;
        }
    }
    return svs;
}

vector<SV> Insdeller::call_svs(const Cluster& cluster, const string& ref) {
    bam1_t *aln = bam_init1() ;
    string region = cluster.chrom + ":" + to_string(cluster.s) + "-" + to_string(cluster.e + 1) ;
    hts_itr_t *itr = sam_itr_querys(read_bamindex[omp_get_thread_num()], read_bamhdr[omp_get_thread_num()], region.c_str()) ;
    uint c = 0 ;
    while (sam_itr_next(read_bam[omp_get_thread_num()], itr, aln) > 0) {
        if (aln->core.flag & BAM_FUNMAP || aln->core.flag & BAM_FSUPPLEMENTARY || aln->core.flag & BAM_FSECONDARY) {
            continue ;
        }
        c++ ;
    }
    Fragment rep ;
    for (const Fragment &f : cluster) {
        //cout << f.ref_s << "-" << f.ref_e << endl ;
        if (f.size() > rep.size()) {
            rep = f;
        }
    }
    vector<SV> svs ;
    // find deletions
    // I HATE this
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
                if (r != -1 && last_r != -1) {
                    int svlen = r - last_r - 1 ;
                    if (svlen < 25) {
                        continue ;
                    }
                    string ref_allele = ref.substr(last_r, svlen) ;
                    string alt_allele = ref_allele.substr(0, 1) ;
                    SV sv = SV("DEL", cluster.chrom, last_r, ref_allele, alt_allele, cluster.size(), c, -1, -1) ;
                    cout << sv << endl ;
                    svs.push_back(sv) ; 
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
                if (q != -1 && last_q != -1) {
                    int svlen = q - last_q - 1 ;
                    if (svlen < 25) {
                        continue ;
                    }
                    string ref_allele = ref.substr(last_r, 1) ;
                    string alt_allele = ref_allele + rep.seq.substr(last_q + 1, svlen) ;
                    SV sv = SV("INS", cluster.chrom, last_r, ref_allele, alt_allele, cluster.size(), c, -1, -1) ;
                    cout << sv << endl ;
                    svs.push_back(sv) ; 
                }
            }
        }
        last_q = q ;
        last_r = r ; 
    }
    rep.aligned_pairs.clear() ;
    free(aln) ;
    return svs ;
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

// Old python code:
// 1. Load all SFS and cluster based on type (tcluster)
// 2. Merge close SFS on each read
// 3. Cluster based on position (pcluster)

void Insdeller::call(const string &chrom_seq, ofstream &osam) {
    // algorithm
    // 1. Put fragments that are close to each other togethe
    // 2. Separate fragments based on type, all DELs go together, all INS go together, fragments that have both go together
    // TODO: does this make sense? why not just put all close fragments together?
    // 3. (why) merge close fragments withe each other, how close? what happens if there is an INS between two DELs
    // 4. Extend each fragment to reach a unique 6bp sequence
    // 5. Cluster based on 6bp sequences, if we have one cluster, then easy, just genotype it
    //      a. if we have one cluster then easy just genotype it
    //      b. if we have multiple clusters, cluster them based on length and overlap
    //          1. if we end up with one or two clusters then predicts SVs
    //          2. if we end up with more, then pass to POA or miniasm

    vector<Cluster> position_clusters = position_cluster() ;
    cout << "Sorted reads into " << position_clusters.size() << " P-clusters." << endl ;

    vector<vector<SV>> _svs(config->threads) ;

    int n = 0 ;
    vector<int> cluster_status[config->threads] ;
    for (int i = 0; i < config->threads; i++) {
        cluster_status[i].reserve(100 + 1) ;
    }

    int m = 0 ;
    #pragma omp parallel for num_threads(config->threads) schedule(dynamic, 1)
    for (int i = 0; i < position_clusters.size(); i++) {
        auto& pc = position_clusters[i] ;
        if (pc.size() < 2) {
            continue ;
        }

        //Parsoa: I don't like this
        vector<Cluster> type_clusters = type_cluster(pc);
        for (Cluster &tc : type_clusters) {
            if (tc.size() < 2) {
                continue ;
            }
            //if (tc.s != 37667602) {
            //    continue ;
            //}
            //if (extended_tc.size() < 2) {
            //    continue ;
            //}
            vector<SV> svs ;
            vector<SV> cl_svs ;
            auto breakpoints = cluster_breakpoints(tc, 0.7) ;
            if (breakpoints.size() == 1 && tc.get_type() != 2) {
                cout << "Using normal genotyper" << endl ;
                svs = call_svs(breakpoints[0], chrom_seq) ;
            } else {
                // extend and pass to POA
                for (auto breakpoint: breakpoints) {
                    auto __svs = call_poa_svs(breakpoint, chrom_seq, osam) ;
                    svs.insert(svs.begin(), __svs.begin(), __svs.end()) ;
                }
            }
            //if (config->assemble) {
            //    svs = haplotyper.haplotype(extended_tc) ;
            //} else {
            //    int n = cluster_anchors(breakpoint) ; 
            //}
            //for (auto& breakpoint: breakpoints) {
            //    if (n == 1) {
            //        // we can call SV
            //        svs = call_svs(breakpoint) ;
            //        cout << "Heterzoygous" << endl ;
            //    } else {
            //        auto sclusters = scluster(breakpoint) ;
            //        if (sclusters.size() == 1) {
            //            svs = call_svs(breakpoint) ;
            //            cout << "Heterzoygous" << endl ;
            //        } else {
            //            lprint({"Delegating to POA.."}) ;
            //            svs = haplotyper.haplotype(extended_tc) ;
            //        }
            //    }
            //}
            cl_svs.insert(cl_svs.end(), svs.begin(), svs.end()) ;
            vector<SV> usvs = remove_duplicate_svs(cl_svs) ;
            for (const SV &sv : usvs) {
                if (abs(sv.l) >= 25 && (sv.ngaps <= 2 || (sv.ngaps > 2 && sv.w > 10))) {
                    _svs[omp_get_thread_num()].push_back(sv);
                }
            }
        }
        if (i % 100 == 0) {
            cout << "=================================================== " << endl ;
            cout << "[P] " << i << " completed.." << endl ;
        }
    }
    cout << "Done.." << endl ;
    cout << n << endl ;
    //for (int i = 0; i < 100; i++) {
    //    int m = 0 ;
    //    for (int j = 0; j < config->threads; j++) {
    //        m += cluster_status[j][i] ;
    //    }
    //    cout << i << " : " << m << endl ;
    //}
    vector<SV> merged_svs ; 
    for (int i = 0; i < _svs.size(); i++) {
        for (const SV &sv : _svs[i]) {
            merged_svs.push_back(sv);
            //vartree.insert({sv.s - 1000, sv.e + 1000});
        }
    }
    osvs = remove_duplicate_svs(merged_svs) ;
    cout << "Discovered " << osvs.size() << " SVs.. " << endl ;
    
    //vartree.deoverlap();
    for (int i = 0; i < config->threads; i++) {
        sam_close(read_bam[i]) ;
    }
}
