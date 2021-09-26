#include "extender.hpp"

using namespace std ;
using namespace lib_interval_tree;

Extender::Extender(const string& _chrom, unordered_map<string, vector<SFS>>* _SFSs) {
    SFSs = _SFSs ;
    chrom = _chrom ;
    ref_seq = string(chromosome_seqs[chrom]) ;
    config = Configuration::getInstance() ;
    bam_file = hts_open(config->bam.c_str(), "r") ;
    bam_index = sam_index_load(bam_file, config->bam.c_str()) ;
    bam_header = sam_hdr_read(bam_file) ;
}

pair<int, int> Extender::get_unique_kmers(const vector<pair<int, int>> &alpairs, const uint k, const bool from_end) {
    if (alpairs.size() < k) {
        return make_pair(-1, -1) ;
    }

    map<string, int> kmers ;
    string kmer_seq ;
    int i = 0 ;
    while (i < alpairs.size() - k + 1) {
        bool skip = false;
        for (int j = i; j < i + k; j++) {
            if (alpairs[j].first == -1 || alpairs[j].second == -1) {
                // we want clean kmers only - ie placed kmers, no insertions or deletions
                skip = true ;
                i = j + 1 ; // jump to next possible start position
                break ;
            }
        }
        if (skip) {
            continue ;
        }
        kmer_seq = ref_seq.substr(alpairs[i].second, k);
        // note: after read reconstruction, if the kmer is placed, than the kmer on the read and the kmer on the refernece should be the same
        kmers[kmer_seq] += 1;
    }
    pair<int, int> last_kmer = make_pair(-1, -1) ;
    i = 0 ;
    while (i < alpairs.size() - k + 1) {
        int offset = i ;
        if (from_end) {
            offset = alpairs.size() - k - i ;
        }
        bool skip = false;
        for (int j = offset; j < offset + k; j++) {
            if (alpairs[j].first == -1 || alpairs[j].second == -1) {
                skip = true ;
                i = j + 1 ; // jump to next possible start position
                break ;
            }
        }
        if (skip) {
            continue ;
        }
        last_kmer = alpairs[offset] ;
        kmer_seq = ref_seq.substr(alpairs[offset].second, k);
        if (kmers[kmer_seq] == 1) {
            break ;
        }
    }
    return last_kmer ;
}

void Extender::extend() {
    bam1_t *aln = bam_init1() ;
    hts_itr_t *itr = sam_itr_querys(bam_index, bam_header, chrom.c_str()) ;
    uint n_al = 0 ;
    while (sam_itr_next(bam_file, itr, aln) > 0) {
        if (aln->core.flag & BAM_FUNMAP || aln->core.flag & BAM_FSUPPLEMENTARY || aln->core.flag & BAM_FSECONDARY) {
            continue ;
        }

        char *qname = bam_get_qname(aln);
        if (SFSs->find(qname) == SFSs->end()) {
            continue;
        }

        uint32_t *cigar = bam_get_cigar(aln);

        vector<pair<int, int>> alpairs = get_aligned_pairs(aln);

        // NOTE: we may have more sfs on a clipped read, but all of them will produce the same clip
        pair<uint, uint> lclip ;
        pair<uint, uint> rclip ;
        int last_pos = 0 ;
        for (SFS &sfs : SFSs->at(qname)) {
            int s = sfs.s;
            int e = sfs.s + sfs.l - 1;
            int aln_start = -1 ;
            int aln_end = -1 ;
            vector<pair<int, int>> local_alpairs;
            // find start and end of SFS ibn read's alignment
            int refs = -1 ;
            int refe = -1 ;
            for (int i = last_pos; i < alpairs.size(); i++) {
                int q = alpairs[i].first;
                int r = alpairs[i].second;
                if (q == -1 || r == -1) {
                    continue ;
                }
                else if (q < s) { // <= seems more correct to me but using < we are more flexible
                    last_pos = i ;
                    refs = r ;
                    aln_start = i ;
                }
                else if (q > e) { // >= seems more correct but using > we are more flexible
                    refe = r ;
                    aln_end = i ;
                    break ;
                }
            }
            // We extract the local alignment of the region of interest
            if (refs == -1 && refe == -1) {
                // we couldn't place the first and the last base, so we skip this - otherwise we'll end up considering the entire read
                ++skip_1 ;
                continue ;
            } else if (refs == -1) {
                uint op = bam_cigar_op(*(cigar + 0));
                uint l = bam_cigar_oplen(*(cigar + 0));
                if (op == BAM_CSOFT_CLIP) {
                    lclip = make_pair(aln->core.pos, l);
                }
                // we skip this SFS
                continue ;
            } else if (refe == -1) {
                uint op = bam_cigar_op(*(cigar + aln->core.n_cigar - 1));
                uint l = bam_cigar_oplen(*(cigar + aln->core.n_cigar - 1));
                if (op == BAM_CSOFT_CLIP) {
                    rclip = make_pair(bam_endpos(aln), l);
                }
                // we skip this SFS
                continue ;
            } else {
                // we placed the first and last base, so we extract the alignment (ie a substring)
                int last_r = refs - 1;
                for (uint i = aln_start; i <= aln_end; ++i) {
                    int q = alpairs[i].first ;
                    int r = alpairs[i].second ;
                    if (r == -1) {
                        if (refs <= last_r && last_r <= refe) {
                            local_alpairs.push_back(make_pair(q, r));
                        }
                    } else {
                        last_r = r;
                        if (refs <= r && r <= refe) {
                            local_alpairs.push_back(make_pair(q, r));
                        }
                    }
                    // We break when we found a placed base at or after the reference end
                    if (q != -1 && r != -1 && r >= refe) {
                        break ;
                    }
                }
            }

            // only if we have been able to place the SFS..
            // 3 ..we extract the maxw (100) pairs preceding the region of interest
            vector<pair<int, int>> pre_alpairs ;
            uint n = 0 ;
            for (int i = aln_start - 1; i >= 0; --i) {
                int q = alpairs[i].first;
                int r = alpairs[i].second;
                if (n < maxw) {
                    pre_alpairs.push_back(make_pair(q, r));
                    n++ ;
                }
            }
            reverse(pre_alpairs.begin(), pre_alpairs.end());
            // 4 we extract the maxw (100) pairs following the region of interest
            vector<pair<int, int>> post_alpairs;
            n = 0;
            for (uint i = aln_end + 1; i < alpairs.size(); ++i) {
                int q = alpairs[i].first;
                int r = alpairs[i].second;
                if (n < maxw) {
                    post_alpairs.push_back(make_pair(q, r));
                    n++ ;
                }
            }

            // 5 we get the unique kmer in the upstream and downstream maxw-bp regions
            pair<int, int> prekmer = get_unique_kmers(pre_alpairs, kmer_size, true);    // true for first kmer found (shorter cluster)
            pair<int, int> postkmer = get_unique_kmers(post_alpairs, kmer_size, false); // false for first kmer found (shorter cluster)

            // if we couldn't place a kmer, we just get the entire region
            if (prekmer.first == -1 || prekmer.second == -1) {
                prekmer.first = local_alpairs.front().first ;
                prekmer.second = local_alpairs.front().second ;
            }
            if (postkmer.first == -1 || postkmer.second == -1) {
                postkmer.first = local_alpairs.back().first ;
                postkmer.second = local_alpairs.back().second ;
            }
            // if also the entire region is not correctly placed, then we skip it
            // NOTE: I think we can solve this by increasing maxw
            if (prekmer.first == -1 || prekmer.second == -1 || postkmer.first == -1 || postkmer.second == -1) {
                ++skip_2 ;
                continue ;
            }
            // FIXME: understand why this is happening (chr16 on full giab genome)
            if ((uint) prekmer.second > postkmer.second + kmer_size) {
                cerr << "Error on " << qname << ". SFS starting at " << sfs.s << " (length " << sfs.l << ")." << endl;
            } else {
                extended_sfs.push_back(ExtSFS(string(chrom), string(qname), prekmer.second, postkmer.second + kmer_size));
            }
            if (lclip.second > 0) {
                clips.push_back(Clip(qname, lclip.first, lclip.second, true)) ;
            }
            if (rclip.second > 0) {
                clips.push_back(Clip(qname, rclip.first, rclip.second, false)) ;
            }
            ++n_al ;
            if (n_al % 10000 == 0) {
                cerr << "Parsed " << n_al << " alignments." << endl;
            }
        }
    }
}

void Extender::cluster() {
    interval_tree_t<int> tree;
    for (const ExtSFS &sfs : extended_sfs) {
        vector<pair<int, int>> overlaps ;
        tree.overlap_find_all({sfs.s, sfs.e}, [&overlaps](auto iter) {
            overlaps.push_back(make_pair(iter->low(), iter->high()));
            return true;
        });
        if (overlaps.empty()) {
            tree.insert({sfs.s, sfs.e});
        } else {
            int mins = sfs.s;
            int maxe = sfs.e;
            for (const pair<int, int> overlap : overlaps) {
                mins = min(mins, overlap.first);
                maxe = max(maxe, overlap.second);
            }
            tree.insert({mins, maxe});
        }
    }
    tree.deoverlap();

    map<pair<int, int>, vector<ExtSFS>> _clusters ;
    for (const ExtSFS &sfs : extended_sfs) {
        auto overlap = tree.overlap_find({sfs.s, sfs.e});
        _clusters[make_pair(overlap->low(), overlap->high())].push_back(sfs);
    }

    // --- CLUSTERS CLEANING
    cerr << "Analyzing " << _clusters.size() << " clusters from " << extended_sfs.size() << " extSFSs.." << endl;

    for (const auto &cluster : _clusters) {
        set<string> reads;
        int cluster_s = numeric_limits<int>::max() ;
        int cluster_e = 0;
        for (const ExtSFS &esfs: cluster.second) {
            cluster_s = min(cluster_s, esfs.s);
            cluster_e = max(cluster_e, esfs.e);
            reads.insert(esfs.qname);
        }
        uint cluster_size = reads.size();
        if (cluster_size < minsupp) {
            ++small_cl ;
            continue ;
        }

        Cluster global_cluster = Cluster(chrom, cluster_s, cluster_e);

        uint cov = 0;
        bam1_t* aln = bam_init1();
        string region = chrom + ":" + to_string(cluster_s) + "-" + to_string(cluster_e);
        hts_itr_t *itr = sam_itr_querys(bam_index, bam_header, region.c_str());
        char* seq = (char*) malloc(10000) ;
        uint32_t len = 0 ;
        // add second level of parallelization inside Extender
        while (sam_itr_next(bam_file, itr, aln) > 0) {
            if (aln->core.flag & BAM_FUNMAP || aln->core.flag & BAM_FSUPPLEMENTARY || aln->core.flag & BAM_FSECONDARY) {
                continue;
            }
            ++cov; // FIXME: this cov takes into account also reads starting or ending inside the cluster (maybe we should skip those?)

            char *qname = bam_get_qname(aln);
            if (reads.find(qname) == reads.end()) {
                continue;
            }

            uint32_t l = aln->core.l_qseq;
            if (l > len) {
                if (len > 0) {
                    free(seq) ;
                }
                len = l ;
                seq = (char*) malloc(l + 1) ;
            }
            uint8_t *q = bam_get_seq(aln) ; //quality string
            for (int i = 0; i < l; i++){
                seq[i] = seq_nt16_str[bam_seqi(q, i)]; //gets nucleotide id and converts them into IUPAC id.
            }
            seq[l] = '\0';

            vector<pair<int, int>> alpairs = get_aligned_pairs(aln);
            int qs = -1, qe = -1;
            // getting starting and ending positions on read sequence aligning to start/end of cluster
            for (int i = alpairs.size() - 1; i >= 0; --i) {
                // finding starting position
                int q = alpairs[i].first;
                int r = alpairs[i].second;
                if (q == -1 || r == -1) {
                    continue;
                }
                if (r <= cluster_s) {
                    qs = q;
                    break;
                }
            }
            for (uint i = 0; i < alpairs.size(); ++i) {
                // finding ending position
                int q = alpairs[i].first;
                int r = alpairs[i].second;
                if (q == -1 || r == -1) {
                    continue;
                }
                if (r >= cluster_e) {
                    qe = q;
                    break;
                }
            }
            if (qs == -1 || qe == -1) {
                // reads starts or ends inside the cluster
                // TODO: get only remaining prefix/suffix? but this may make POA and realignment harder
                ++skip_3;
            } else {
                string seq(seq, qs, qe - qs + 1);
                global_cluster.add(seq);
            }
        }
        if (global_cluster.size() >= minsupp) {
            global_cluster.set_cov(cov) ;
            clusters.push_back(global_cluster) ;
            ++extcl;
        } else {
            ++small_extcl;
        }
    }
    cerr << "Dumped " << extcl << " clusters." << endl << endl;
    cerr << "Skipped due to unplaced alignment: " << skip_1 << endl;
    cerr << "Skipped due to unplaced extension: " << skip_2 << endl;
    cerr << "Skipped due to start/end inside extended cluster: " << skip_3 << endl;
    cerr << "Skipped extended clusters (due to support): " << small_extcl << endl;
}

vector<Cluster> Extender::cluster_by_length(const Cluster& cluster) {
    vector<Cluster> clusters_by_len;
    for (const string &seq : cluster.get_seqs()) {
        uint i;
        for (i = 0; i < clusters_by_len.size(); ++i) {
            if (abs((int)clusters_by_len[i].get_len() - (int)seq.size()) <= mind) {
                break;
            }
        }
        if (i == clusters_by_len.size()) {
            clusters_by_len.push_back(Cluster(cluster.chrom, cluster.s, cluster.e, cluster.cov));
        }
        clusters_by_len[i].add(seq);
    }
    return clusters_by_len ;
}

vector<pair<uint, char>> Extender::parse_cigar(string cigar) {
    // -- Parsing CIGAR
    vector<pair<uint, char>> cigar_pairs ;
    regex r("([0-9]+)([MIDNSHPX=])") ;
    regex_iterator<string::iterator> rit(cigar.begin(), cigar.end(), r) ;
    regex_iterator<string::iterator> rend ;
    uint nv = 0 ;
    while (rit != rend) {
        int l = stoi(rit->str().substr(0, rit->str().size() - 1)) ;
        char op = rit->str().substr(rit->str().size() - 1, 1)[0] ;
        cigar_pairs.push_back(make_pair(l, op)) ;
        if (l > 25 && (op == 'I' || op == 'D')) {
            ++nv ;
        }
        ++rit ;
    }
    return cigar_pairs ;
}

void Extender::call() {
    for (const Cluster &cluster: clusters) {
        if (cluster.size() < minw) {
            continue;
        }
        // --- Splitting cluster by sequence length
        const auto& clusters_by_len = cluster_by_length(cluster) ;
        // --- Sorting clusters by #sequences to get first 2 most weighted clusters
        int i_max1 = -1 ;
        int i_max2 = -1 ;
        uint v_max1 = 0 ;
        uint v_max2 = 0 ;
        for (uint i = 0; i < clusters_by_len.size(); ++i) {
            if (clusters_by_len[i].size() > v_max1) {
                v_max2 = v_max1;
                i_max2 = i_max1;
                v_max1 = clusters_by_len[i].size();
                i_max1 = i;
            }
            else if (clusters_by_len[i].size() > v_max2) {
                v_max2 = clusters_by_len[i].size();
                i_max2 = i;
            }
        }
        vector<int> maxs ({i_max1, i_max2});
        // Genotyping the two most-weighted clusters
        for (const int i : maxs) {
            if (i == -1) {
                continue;
            }
            Cluster c = clusters_by_len[i];
            if (c.size() < minw) {
                continue;
            }
            // --- Local realignment
            vector<SV> svs; // svs on current cluster
            string ref = ref_seq.substr(c.s, c.e - c.s + 1);
            string consensus = c.poa() ;
            parasail_result_t *result = NULL;
            result = parasail_nw_trace_striped_16(consensus.c_str(), consensus.size(), ref.c_str(), ref.size(), 10, 1, &parasail_blosum62);
            parasail_cigar_t *cigar = parasail_result_get_cigar(result, consensus.c_str(), consensus.size(), ref.c_str(), ref.size(), NULL);
            string cigar_str = parasail_cigar_decode(cigar);
            int score = result->score ;
            parasail_cigar_free(cigar) ;
            parasail_result_free(result) ;
            alignments.push_back(Consensus(consensus, cigar_str, chrom, c.s, c.e)) ;
            // -- Extracting SVs
            uint rpos = c.s ; // position on reference
            uint cpos = 0 ;   // position on consensus
            auto cigar_pairs = parse_cigar(cigar_str) ; 
            for (const auto cigar_pair: cigar_pairs) {
                uint l = cigar_pair.first;
                char op = cigar_pair.second;
                if (op == '=' || op == 'M') {
                    rpos += l;
                    cpos += l;
                } else if (op == 'I') {
                    if (l > 25) {
                        svs.push_back(SV("INS", c.chrom, rpos, ref_seq.substr(rpos - 1, 1), consensus.substr(cpos, l), c.size(), c.cov, l, 0, score));
                    }
                    cpos += l;
                } else if (op == 'D') {
                    if (l > 25) {
                        svs.push_back(SV("DEL", c.chrom, rpos, ref_seq.substr(rpos - 1, l), ref_seq.substr(rpos - 1, 1), c.size(), c.cov, -l, 0, score));
                    }
                    rpos += l;
                }
            }
            // --- combine SVs on same consensus ---
            vector<SV> merged_svs;
            SV ins;
            SV del;
            for (const auto &sv : svs) {
                if (sv.l < 0) {
                    if (del.l == 0) {
                        del = sv ;
                        continue ;
                    } else {
                        // TODO: Parsoa: isn't this wrong? what if they the deletions are not overlapping or what if there is an ins in between?
                        del.l += sv.l ;
                        del.refall += sv.refall ;
                        del.e += (-sv.l) ;
                    }
                } else {
                    if (ins.l == 0) {
                        ins = sv;
                        continue;
                    } else {
                        ins.l += sv.l;
                        ins.altall += sv.altall;
                    }
                }
            }
            if (ins.l != 0) {
                merged_svs.push_back(ins);
            }
            if (del.l != 0) {
                merged_svs.push_back(del);
            }
            // -- Combine svs with same length (maybe useless now - only if diploid mode)
            vector<SV> comb_svs;
            for (const SV &msv : merged_svs) {
                bool newsv_flag = true;
                for (SV &csv : comb_svs) {
                    if (abs(csv.l - msv.l) <= 10) {
                        csv.w += msv.w;
                        newsv_flag = false;
                    }
                }
                if (newsv_flag) {
                    comb_svs.push_back(msv);
                }
            }
            for (const SV &sv : comb_svs) {
                this->svs.push_back(sv);
            }
        }
    }
}
