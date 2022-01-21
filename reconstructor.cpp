#include "reconstructor.hpp"

using namespace std ;

//typedef struct bam1_t {
//    bam1_core_t core; // won't change
//    uint64_t id;      // won't change
//    uint8_t *data;    // has to reallocate
//    int l_data;       // will change
//    uint32_t m_data;  // won't change
//    uint32_t mempolicy:2, :30 /* Reserved */;
//} bam1_t;

//typedef struct bam1_core_t {
//    hts_pos_t pos;      // won't change 
//    int32_t tid;        // won't change
//    uint16_t bin;       // just copy 
//    uint8_t qual;       // won't change
//    uint8_t l_extranul; // won't change 
//    uint16_t flag;      // won't change
//    uint16_t l_qname;   // won't change
//    uint32_t n_cigar;   // will change
//    int32_t l_qseq;     // will change
//    int32_t mtid;       // just copy 
//    hts_pos_t mpos;     // just copy 
//    hts_pos_t isize;    // may or may ont change?
//} bam1_core_t;

int do_realloc_bam_data(bam1_t *b, size_t desired) {
    //cout << "Extending BAM entry." << endl ;
    uint32_t new_m_data ;
    uint8_t *new_data ;
    new_m_data = desired ;
    kroundup32(new_m_data) ;
    if (new_m_data < desired) {
        errno = ENOMEM ; // Not strictly true but we can't store the size
        return -1 ;
    }
    new_data = (uint8_t*) realloc(b->data, new_m_data) ;
    if (!new_data) return -1 ;
    b->data = new_data ;
    b->m_data = new_m_data ;
    return 0 ;
}

int realloc_bam_data(bam1_t *b, size_t desired) {
    if (desired <= b->m_data) return 0 ;
    return do_realloc_bam_data(b, desired) ;
}

void rebuild_bam_entry(bam1_t* alignment, char* seq, uint8_t* qual, vector<pair<uint32_t, uint32_t>> cigar) {
    auto l_aux = bam_get_l_aux(alignment) ;
    uint8_t* aux = (uint8_t*) malloc(sizeof(uint8_t) * l_aux) ;
    memcpy(aux, alignment->data + alignment->l_data - l_aux, l_aux) ;
    // update core
    alignment->core.n_cigar = cigar.size() ;
    alignment->core.l_qseq = strlen(seq) ;
    int l = strlen(seq) ;
    // rebuild data
    int l_data = alignment->core.l_qname + (4 * alignment->core.n_cigar) + ((l + 1) >> 1) + l + l_aux ;
    realloc_bam_data(alignment, l_data) ;
    alignment->l_data = l_data ;
    // copy qname
    int offset = alignment->core.l_qname ;
    // copy cigar
    uint8_t* cigar_encoded = encode_cigar(cigar) ;
    memcpy(alignment->data + offset, cigar_encoded, 4 * alignment->core.n_cigar) ;
    offset += 4 * alignment->core.n_cigar ;
    free(cigar_encoded) ;
    // copy seq data - have to convert seq
    uint8_t* seq_bytes = encode_bam_seq(seq) ;
    memcpy(alignment->data + offset, seq_bytes, (l + 1) >> 1) ;
    free(seq_bytes) ;
    offset += ((l + 1) >> 1) ;
    // copy quality
    memcpy(alignment->data + offset, qual, l) ;
    offset += l ;
    // copy aux
    memcpy(alignment->data + offset, aux, l_aux) ;
    free(aux) ;
}

void Reconstructor::reconstruct_read(bam1_t* alignment, char* read_seq, string chrom, int _i, int _j, int _k) {
    auto cigar_offsets = decode_cigar(alignment) ;
    int l = 0 ;
    // try and filter unintenresting reads early on
    bool should_ignore = true ;
    for (auto op: cigar_offsets) {
        l += op.first ;
    }
    if (read_seq_max_lengths[_i][_j][_k] < l) {
        free(read_seqs[_i][_j][_k]) ;
        free(new_read_seqs[_i][_j][_k]) ;
        free(new_read_quals[_i][_j][_k]) ;
        //
        read_seqs[_i][_j][_k] = (char*) malloc(sizeof(char) * (l + 1)) ;
        new_read_seqs[_i][_j][_k] = (char*) malloc(sizeof(char) * (l + 1)) ;
        new_read_quals[_i][_j][_k] = (uint8_t*) malloc(sizeof(char) * (l + 1)) ;
        read_seq_max_lengths[_i][_j][_k] = l ;
    }
    //
    int n = 0 ;
    int m = 0 ;
    int ref_offset = alignment->core.pos ;
    int ins_offset = 0 ;
    int del_offset = 0 ;
    int match_offset = 0 ;
    int soft_clip_offset = 0 ;
    char* new_seq = new_read_seqs[_i][_j][_k] ; 
    uint8_t* qual = bam_get_qual(alignment) ;
    uint8_t* new_qual = new_read_quals[_i][_j][_k] ; 
    int pos = alignment->core.pos + 1 ; // this is 0-based, variant cpoordinates are 1-based
    // Modify current bam1_t* struct
    auto& core = alignment->core ;
    vector<pair<uint32_t, uint32_t>> new_cigar ;
    int m_diff = 0 ;
    double num_match = 0 ;
    double num_mismatch = 0 ;
    char* ref_seq = chromosome_seqs[chrom] ;
    while (true) {
        if (m == cigar_offsets.size()) {
            break ;
        }
        if (cigar_offsets[m].second == BAM_CMATCH || cigar_offsets[m].second == BAM_CEQUAL || cigar_offsets[m].second == BAM_CDIFF) {
            memcpy(new_seq + n, ref_seq + ref_offset, cigar_offsets[m].first) ;
            memcpy(new_qual + n, qual + match_offset + ins_offset + soft_clip_offset, cigar_offsets[m].first) ;
            n += cigar_offsets[m].first ;
            for (int j = 0; j < cigar_offsets[m].first; j++) {
                num_mismatch += 1 ? ref_seq[ref_offset + j] != read_seq[match_offset + ins_offset + soft_clip_offset + j] : 0 ;
            }
            num_match += cigar_offsets[m].first ;
            ref_offset += cigar_offsets[m].first ;
            match_offset += cigar_offsets[m].first ;
            if (new_cigar.size() >= 1 && new_cigar[new_cigar.size() - 1].second == BAM_CMATCH) {
                new_cigar[new_cigar.size() - 1].first += cigar_offsets[m].first + m_diff ;
            } else {
                new_cigar.push_back(make_pair(cigar_offsets[m].first + m_diff, BAM_CMATCH)) ;
            }
            m_diff = 0 ;
        } else if (cigar_offsets[m].second == BAM_CINS) {
            if (cigar_offsets[m].first <= config->min_indel_length) {
                // if a short INDEL then just don't add it to read
            } else {
                // for long INS, this is probably a SV so add it to the read
                should_ignore = false ;
                memcpy(new_seq + n, read_seq + soft_clip_offset + match_offset + ins_offset, cigar_offsets[m].first) ;
                memcpy(new_qual + n, qual + soft_clip_offset + match_offset + ins_offset, cigar_offsets[m].first) ;
                n += cigar_offsets[m].first ;
                new_cigar.push_back(cigar_offsets[m]) ;
            }
            ins_offset += cigar_offsets[m].first ;
        } else if (cigar_offsets[m].second == BAM_CDEL) {
            if (cigar_offsets[m].first <= config->min_indel_length) {
                // if a short DEL so let's just fix it
                memcpy(new_seq + n, ref_seq + ref_offset, cigar_offsets[m].first) ;
                memcpy(new_qual + n, qual + soft_clip_offset + match_offset + ins_offset, cigar_offsets[m].first) ;
                n += cigar_offsets[m].first ;
                m_diff += cigar_offsets[m].first ;
            } else {
                // for long DEL, this is probably a SV so let it be what it was
                should_ignore = false ;
                new_cigar.push_back(cigar_offsets[m]) ;
            }
            del_offset += cigar_offsets[m].first ;
            ref_offset += cigar_offsets[m].first ;
        } else if (cigar_offsets[m].second == BAM_CSOFT_CLIP) {
            should_ignore = false ;
            memcpy(new_seq + n, read_seq + soft_clip_offset + match_offset + ins_offset, cigar_offsets[m].first) ;
            memcpy(new_qual + n, qual + soft_clip_offset + match_offset + ins_offset, cigar_offsets[m].first) ;
            n += cigar_offsets[m].first ;
            soft_clip_offset += cigar_offsets[m].first ;
            new_cigar.push_back(cigar_offsets[m]) ;
        } else {
            cout << "Illegal Cigar OP" << endl ;
            break ;
            //if (cigar_offsets[m].second == BAM_CPAD || cigar_offsets[m].second == BAM_CHARD_CLIP || cigar_offsets[m].second == BAM_CBACK) {
        }
        m += 1 ;
    }
    new_seq[n] = '\0' ;
    new_qual[n] = '\0' ;
    //int r = 0 ;
    //for (auto op: new_cigar) {
    //    if (op.second != BAM_CDEL) {
    //        n -= op.first ;
    //    }
    //    r += op.first ;
    //}
    //assert(n == 0) ;
    char *qname = bam_get_qname(alignment) ;
    //cout << n << " " << strlen(new_seq) << " " << r << endl ;
    //cout << bam_get_qname(alignment) << endl ;
    // only do this on first processing thread
    if (omp_get_thread_num() == 2) {
        global_num_bases += num_match ;
        global_num_mismatch += num_mismatch ;
    }
    // how many errors and SNPs do we expect? 1/1000 each, so say if we see more than twice that then don't correct
    if (config->selective) {
        if (num_mismatch / num_match > 3 * expected_mismatch_rate) {
            if (omp_get_thread_num() == 3) {
                num_ignored_reads += 1 ;
            }
            return ;
        }
        // if we have so many deletions and insertions, then abort
        if (ins_offset + del_offset > 0.7 * strlen(read_seq)) {
            return ;
        }
    }
    if (should_ignore) {
        // check mismatch rate
        ignored_reads[omp_get_thread_num() - 2].push_back(qname) ;
        return ;
    }
    reconstructed_reads[omp_get_thread_num() - 2].push_back(qname) ;
    rebuild_bam_entry(alignment, new_seq, new_qual, new_cigar) ;
}

void Reconstructor::process_batch(vector<bam1_t*> bam_entries, int p, int i) {
    bam1_t* alignment ;
    for (int b = 0; b < bam_entries.size(); b++) {
        alignment = bam_entries[b] ;
        if (alignment == nullptr) {
            break ;
        }
        if (alignment->core.flag & BAM_FUNMAP || alignment->core.flag & BAM_FSUPPLEMENTARY || alignment->core.flag & BAM_FSECONDARY) {
            continue ;
        }
        if (alignment->core.l_qseq < 2) {
            //cerr << "Read too short, ignoring.." << endl ;
            continue ;
        }
        if (alignment->core.tid < 0) {
            continue ;
        }
        string chrom(bam_header->target_name[alignment->core.tid]) ;
        if (chromosome_seqs.find(chrom) == chromosome_seqs.end()) {
            continue ;
        }
        reconstruct_read(alignment, read_seqs[p][i][b], chrom, p, i, b) ;
    }
}

// BAM writing based on https://www.biostars.org/p/181580/
void Reconstructor::run() {
    config = Configuration::getInstance() ;
    load_chromosomes(config->reference) ;
    // parse arguments
    bam_file = hts_open(config->bam.c_str(), "r") ;
    bam_index = sam_index_load(bam_file, config->bam.c_str()) ;
    bam_header = sam_hdr_read(bam_file) ; //read header
    bgzf_mt(bam_file->fp.bgzf, 8, 1) ;
    auto out_bam_path = config->workdir + (config->selective ? "/reconstructed.selective.bam" : "/reconstructed.bam") ;
    out_bam_file = hts_open(out_bam_path.c_str(), "wb") ;
    bgzf_mt(out_bam_file->fp.bgzf, 8, 1) ;
    int r = sam_hdr_write(out_bam_file, bam_header) ;
    if (r < 0) {
        lprint({"Can't write corrected BAM header, aborting.."}, 2);
        return ;
    }
    // allocate stuff
    int modulo = 3 ;
    int batch_size = (10000 / config->threads) * config->threads ;
    for (int i = 0; i < modulo; i++) {
        read_seqs.push_back(vector<vector<char*>>(config->threads)) ; // current and next output
        new_read_seqs.push_back(vector<vector<char*>>(config->threads)) ; // current and next output
        new_read_quals.push_back(vector<vector<uint8_t*>>(config->threads)) ; // current and next output
        read_seq_lengths.push_back(vector<vector<int>>(config->threads)) ; // current and next output
        read_seq_max_lengths.push_back(vector<vector<int>>(config->threads)) ; // current and next output
        for (int j = 0; j < config->threads; j++) {
            for (int k = 0; k < batch_size / config->threads; k++) {
                read_seqs[i][j].push_back((char*) malloc(sizeof(char) * (30001))) ;
                new_read_seqs[i][j].push_back((char*) malloc(sizeof(char) * (30001))) ;
                new_read_quals[i][j].push_back((uint8_t*) malloc(sizeof(uint8_t) * (30001))) ;
                //
                read_seq_lengths[i][j].push_back(30000) ;
                read_seq_max_lengths[i][j].push_back(30000) ;
            }
        }
    }
    for (int i = 0; i < modulo; i++) {
        bam_entries.push_back(vector<vector<bam1_t*>>(config->threads)) ;
        for (int j = 0; j < config->threads; j++) {
            for (int k = 0; k < batch_size / config->threads; k++) {
                bam_entries[i][j].push_back(bam_init1()) ;
            }
        }
    }
    int b = 0 ;
    lprint({"Loading first batch.."});
    load_batch_bam(config->threads, batch_size, 1) ;
    int p = 1 ;
    ignored_reads.resize(config->threads) ;
    reconstructed_reads.resize(config->threads) ;
    time_t t ;
    time(&t) ;
    bool should_load = true ;
    bool should_process = true ;
    bool should_terminate = false ;
    bool loaded_last_batch = false ;
    int reads_written = 0 ;
    //lprint({"Starting main loop.."}) ;
    while (should_process) {
        //lprint({"Beginning batch", to_string(b + 1)});
        if (!should_load) {
            should_process = false ;
        }
        if (loaded_last_batch) {
            should_load = false ;
        }
        #pragma omp parallel for num_threads(config->threads + 2)
        for(int i = 0; i < config->threads + 2; i++) {
            if (i == 0) {
                if (should_load) {
                    loaded_last_batch = !load_batch_bam(config->threads, batch_size, (p + 1) % modulo) ;
                    //if (loaded_last_batch) {
                    //    lprint({"Last input batch loaded."});
                    //} else {
                    //    lprint({"Loaded."});
                    //}
                }
            } else if (i == 1) {
                if (b >= 1) {
                    int ret = 0 ;
                    for (int k = 0; k < batch_size / config->threads; k++) {
                        for (int j = 0; j < config->threads; j++) {
                            if (bam_entries[(p + 2) % modulo][j][k] != nullptr) {
                                auto alignment = bam_entries[(p + 2) % modulo][j][k] ;
                                ret = sam_write1(out_bam_file, bam_header, bam_entries[(p + 2) % modulo][j][k]) ;
                                reads_written += 1 ;
                                if (ret < 0) {
                                    lprint({"Can't write corrected BAM record, aborting.."}, 2);
                                    should_terminate = true ;
                                    break ;
                                }
                            } else {
                                break ;
                            }
                        }
                    }
                }
            } else {
                if (should_process) {
                    process_batch(bam_entries[p][i - 2], p, i - 2) ;
                }
            }
        }
        if (should_terminate) {
            lprint({"Something went wrong, aborting.."}, 2);
            return ;
        }
        //if (!should_load) {
        //    lprint({"Processed last batch of inputs."});
        //}
        //if (!should_process) {
        //    lprint({"Exiting accelerator loop."});
        //    break ;
        //}
        p += 1 ;
        p %= modulo ;
        b += 1 ;
        time_t s ;
        time(&s) ;
        if (s - t == 0) {
            s += 1 ;
        }
        cerr << "[I] Processed batch " << b << ". Reads so far " << reads_processed << ". Reads per second: " << reads_processed / (s - t) << ". Time: " << s - t << "\n" ;
        cerr << "[I] Processed bases: " << uint64_t(global_num_bases) << ", num mismatch: " << uint64_t(global_num_mismatch) << ", mismatch rate: " << global_num_mismatch / global_num_bases << ", ignored reads: " << num_ignored_reads << "\n" ;
        expected_mismatch_rate = global_num_mismatch / global_num_bases ; 
        cerr << "\x1b[A" ;
        cerr << "\x1b[A" ;
    }
    cerr << endl ;
    cerr << endl ;
    lprint({"Done."});
    sam_close(bam_file) ;
    sam_close(out_bam_file) ;
    dump_reconstructed_read_ids() ;
    lprint({"Loaded", to_string(reads_processed), "reads."});
    lprint({"Wrote", to_string(reads_written), "reads."});
}

void Reconstructor::dump_reconstructed_read_ids() {
    lprint({"Dumping reconstructed read ids.."}) ;
    ofstream qname_file(config->workdir + "/reconstructed_reads.txt") ;
    if (qname_file.is_open()) {
        for (int i = 0; i < config->threads; i++) {
            for (const auto& qname: reconstructed_reads[i]) {
                qname_file << qname << endl ;
            }
        }
    } else {
        lprint({"Error openning reconstructed_reads.txt."}, 2) ;
    }
    qname_file.close() ;
    ofstream ignore_file(config->workdir + "/ignored_reads.txt") ;
    if (ignore_file.is_open()) {
        for (int i = 0; i < config->threads; i++) {
            for (const auto& qname: ignored_reads[i]) {
                ignore_file << qname << endl ;
            }
        }
        ignore_file.close() ;
    } else {
        lprint({"Error openning ignored_read.txt."}, 2) ;
    }
}

bool Reconstructor::load_batch_bam(int threads, int batch_size, int p) {
    int n = 0 ;
    int i = 0 ;
    int m = 0 ;
    while (sam_read1(bam_file, bam_header, bam_entries[p][n % threads][i]) >= 0) {
        auto alignment = bam_entries[p][n % threads][i] ;
        if (alignment == nullptr) {
            break ;
        }
        reads_processed += 1 ;
        uint32_t l = alignment->core.l_qseq ; //length of the read
        if (read_seq_max_lengths[p][n % threads][i] < l) {
            free(read_seqs[p][n % threads][i]) ;
            free(new_read_seqs[p][n % threads][i]) ;
            free(new_read_quals[p][n % threads][i]) ;
            //
            read_seqs[p][n % threads][i] = (char*) malloc(sizeof(char) * (l + 1)) ;
            new_read_seqs[p][n % threads][i] = (char*) malloc(sizeof(char) * (l + 1)) ;
            new_read_quals[p][n % threads][i] = (uint8_t*) malloc(sizeof(char) * (l + 1)) ;
            read_seq_max_lengths[p][n % threads][i] = l ;
            m += 1 ;
        }
        read_seq_lengths[p][n % threads][i] = l ;
        uint8_t *q = bam_get_seq(alignment) ;
        for (int _ = 0; _ < l; _++){
            read_seqs[p][n % threads][i][_] = seq_nt16_str[bam_seqi(q, _)] ;
        }
        read_seqs[p][n % threads][i][l] = '\0' ;
        n += 1 ;
        if (n % threads == 0) {
            i += 1 ;
        }
        if (n == batch_size) {
            break ;
        }
    }
    //cout << m << " reallocations.." << endl ;
    // last batch was incomplete
    if (n != batch_size) {
        for (int j = n % threads; j < threads; j++) {
            //cout << "Terminus at " << j << " " << i << endl ;
            for (int _ = i; _ < batch_size / threads; _++) {
                bam_entries[p][j][_] = nullptr ;
            }
        }
        for (int j = 0; j < n % threads; j++) {
            //cout << "Terminus at " << j << " " << i + 1 << endl ;
            for (int _ = i + 1; _ < batch_size / threads; _++) {
                bam_entries[p][j][_] = nullptr ;
            }
        }
    }
    //lprint({"Loaded", to_string(n), "BAM reads.."});
    return n != 0 ? true : false ;
}
