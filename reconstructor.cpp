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
    cout << "Extending BAM entry." << endl ;
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

void Reconstructor::reconstruct_read(bam1_t* alignment, char* read_seq, string chrom) {
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
    char* new_seq = (char*) malloc(sizeof(char) * (l + 1)) ;
    uint8_t* qual = bam_get_qual(alignment) ;
    uint8_t* new_qual = (uint8_t*) malloc(sizeof(char) * (l + 1)) ;
    int pos = alignment->core.pos + 1 ; // this is 0-based, variant cpoordinates are 1-based
    // Modify current bam1_t* struct
    auto& core = alignment->core ;
    vector<pair<uint32_t, uint32_t>> new_cigar ;
    int m_diff = 0 ;
    double num_match = 0 ;
    double num_mismatch = 0 ;
    while (true) {
        if (m == cigar_offsets.size()) {
            break ;
        }
        if (cigar_offsets[m].second == BAM_CMATCH || cigar_offsets[m].second == BAM_CEQUAL || cigar_offsets[m].second == BAM_CDIFF) {
            for (int j = 0; j < cigar_offsets[m].first; j++) {
                new_seq[n] = chromosome_seqs[chrom][ref_offset + j] ;
                new_qual[n] = qual[soft_clip_offset + match_offset + ins_offset + j] ;
                num_mismatch += 1 ? chromosome_seqs[chrom][ref_offset + j] != read_seq[match_offset + ins_offset + soft_clip_offset + j] : 0 ;
                n++ ;
            }
            ref_offset += cigar_offsets[m].first ;
            match_offset += cigar_offsets[m].first ;
            num_match += cigar_offsets[m].first ;
            if (new_cigar.size() >= 1 && new_cigar[new_cigar.size() - 1].second == BAM_CMATCH) {
                new_cigar[new_cigar.size() - 1].first += cigar_offsets[m].first + m_diff ;
            } else {
                new_cigar.push_back(make_pair(cigar_offsets[m].first + m_diff, BAM_CMATCH)) ;
            }
            m_diff = 0 ;
        } else if (cigar_offsets[m].second == BAM_CINS) {
            if (cigar_offsets[m].first <= 10) {
                // if a short INDEL then just don't add it to read
            } else {
                // for long INS, this is probably a SV so add it to the read
                for (int j = 0; j < cigar_offsets[m].first; j++) {
                    new_seq[n] = read_seq[soft_clip_offset + match_offset + ins_offset + j] ;
                    new_qual[n] = qual[soft_clip_offset + match_offset + ins_offset + j] ; // bases are in read
                    n++ ;
                }
                new_cigar.push_back(cigar_offsets[m]) ;
            }
            ins_offset += cigar_offsets[m].first ;
        } else if (cigar_offsets[m].second == BAM_CDEL) {
            if (cigar_offsets[m].first <= 10) {
                // if a short DEL so let's just fix it
                for (int j = 0; j < cigar_offsets[m].first; j++) {
                    new_seq[n] = chromosome_seqs[chrom][ref_offset + j] ;
                    new_qual[n] = qual[soft_clip_offset + match_offset + ins_offset] ; // just use last observed quality
                    n++ ;
                }
                m_diff += cigar_offsets[m].first ;
            } else {
                // for long DEL, this is probably a SV so let it be what it was
                new_cigar.push_back(cigar_offsets[m]) ;
            }
            del_offset += cigar_offsets[m].first ;
            ref_offset += cigar_offsets[m].first ;
        } else if (cigar_offsets[m].second == BAM_CSOFT_CLIP) {
            for (int j = 0; j < cigar_offsets[m].first; j++) {
                new_seq[n] = read_seq[soft_clip_offset + match_offset + ins_offset + j] ;
                new_qual[n] = qual[soft_clip_offset + match_offset + ins_offset + j] ;
                n++ ;
            }
            soft_clip_offset += cigar_offsets[m].first ;
            new_cigar.push_back(cigar_offsets[m]) ;
        } else if (cigar_offsets[m].second == BAM_CREF_SKIP) {
            // won't happen in DNA alignments
        } else {//if (cigar_offsets[m].second == BAM_CPAD || cigar_offsets[m].second == BAM_CHARD_CLIP || cigar_offsets[m].second == BAM_CBACK) {
            // pass
        }
        m += 1 ;
    }
    new_seq[n] = '\0' ;
    new_qual[n] = '\0' ;
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
    rebuild_bam_entry(alignment, new_seq, new_qual, new_cigar) ;
    free(new_seq) ;
    free(new_qual) ;
}

void Reconstructor::process_batch(vector<bam1_t*> bam_entries) {
    char* seq = (char*) malloc(10000) ;
    uint32_t len = 0 ;
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
        // recover sequence
        uint32_t l = alignment->core.l_qseq ; //length of the read
        if (l > len) {
            if (len > 0) {
                free(seq) ;
            }
            len = l ;
            seq = (char*) malloc(l + 1) ;
        }
        uint8_t *q = bam_get_seq(alignment) ; //quality string
        for (int i = 0; i < l; i++){
            seq[i] = seq_nt16_str[bam_seqi(q, i)]; //gets nucleotide id and converts them into IUPAC id.
        }
        seq[l] = '\0' ; // null terminate
        //cout << bam_get_qname(alignment) << " " << bam_header->target_name[alignment->core.tid] << " " << alignment->core.mpos << endl ;
        reconstruct_read(alignment, seq, chrom) ;
    }
    free(seq) ;
}

// BAM writing based on https://www.biostars.org/p/181580/
void Reconstructor::run() {
    config = Configuration::getInstance() ;
    load_chromosomes(config->reference) ;
    // parse arguments
    bam_file = sam_open(config->bam.c_str(), "r") ;
    bam_header = sam_hdr_read(bam_file) ; //read header
    bgzf_mt(bam_file->fp.bgzf, 8, 1) ;
    auto out_bam_path = config->workdir + (config->selective ? "/reconstructed.selective.bam" : "/reconstructed.bam") ;
    out_bam_file = sam_open(out_bam_path.c_str(), "wb") ;
    int r = sam_hdr_write(out_bam_file, bam_header) ;
    if (r < 0) {
        lprint({"Can't write corrected BAM header, aborting.."}, 2);
    }
    // confidence scores 
    for(int i = 0; i < 2; i++) {
        bam_entries.push_back(vector<vector<bam1_t*>>(config->threads)) ;
    }
    int p = 0 ;
    int b = 0 ;
    int batch_size = (10000 / config->threads) * config->threads ;
    lprint({"Loading first batch.."});
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < config->threads; j++) {
            for (int k = 0; k <= batch_size / config->threads; k++) {
                bam_entries[i][j].push_back(bam_init1()) ;
            }
        }
    }
    load_batch_bam(config->threads, batch_size, p) ;
    // main loop
    time_t t ;
    time(&t) ;
    bool should_load = true ;
    bool should_process = true ;
    bool should_terminate = false ;
    bool loaded_last_batch = false ;
    uint64_t u = 0 ;
    int num_reads = 0 ;
    while (true) {
      lprint({"Beginning batch", to_string(b + 1)});
        for (int i = 0 ; i < config->threads ; i++) {
            u += bam_entries[p][i].size() ;
        }
        if (!should_load) {
            should_process = false ;
        }
        if (loaded_last_batch) {
            should_load = false ;
        }
        #pragma omp parallel for num_threads(config->threads + 2)
        for(int i = 0; i < config->threads + 2; i++) {
            if (i == 0) {
                // write previous batch
                if (b >= 1) {
                    // write BAM
                    int ret = 0 ;
                    for (int k = 0; k < batch_size / config->threads; k++) {
                        for (int j = 0; j < config->threads; j++) {
                            if (bam_entries[(p + 1) % 2][j][k] != nullptr) {
                                auto alignment = bam_entries[(p + 1) % 2][j][k] ;
                                ret = sam_write1(out_bam_file, bam_header, bam_entries[(p + 1) % 2][j][k]);
                                num_reads++ ;
                                if (ret < 0) {
                                    lprint({"Can't write corrected BAM record, aborting.."}, 2);
                                    should_terminate = true ;
                                }
                            } else {
                                break ;
                            }
                        }
                    }
                }
                // load next batch of entries
                if (should_load) {
                    loaded_last_batch = !load_batch_bam(config->threads, batch_size, (p + 1) % 2) ;
                    if (loaded_last_batch) {
                        lprint({"Last input batch loaded."});
                    } else {
                        lprint({"Loaded."});
                    }
                }
            } else if (i == 1) {
                // merge output of previous batch
            } else {
                // process current batch
                if (should_process) {
                    process_batch(bam_entries[p][i - 2]) ;
                }
            }
        }
        if (should_terminate) {
            lprint({"Something went wrong, aborting.."}, 2);
        }
        if (!should_load) {
            lprint({"Processed last batch of inputs."});
        }
        if (!should_process) {
            break ;
        }
        p += 1 ;
        p %= 2 ;
        b += 1 ;
        time_t s ;
        time(&s) ;
        if (s - t == 0) {
            s += 1 ;
        }
        cerr << "[I] Processed batch " << std::left << std::setw(10) << b << ". Reads so far " << std::right << std::setw(12) << u << ". Reads per second: " <<  u / (s - t) << ". Time: " << std::setw(8) << std::fixed << s - t << "\n" ;
        cerr << "[I] Process bases: " << std::left << std::setw(16) << uint64_t(global_num_bases) << ", num mismatch: " << std::setw(16) << uint64_t(global_num_mismatch) << ", mismatch rate: " << global_num_mismatch / global_num_bases << ", ignored reads: " << num_ignored_reads << endl ;
        expected_mismatch_rate = global_num_mismatch / global_num_bases ; 
    }
    lprint({"Done."});
    sam_close(bam_file) ;
    sam_close(out_bam_file) ;
    lprint({"Wrote", to_string(num_reads), "reads."});
}

bool Reconstructor::load_batch_bam(int threads, int batch_size, int p) {
    int n = 0 ;
    int i = 0 ;
    while (sam_read1(bam_file, bam_header, bam_entries[p][n % threads][i]) >= 0) {
        n += 1 ;
        if (n % threads == 0) {
            i += 1 ;
        }
        if (n == batch_size) {
            break ;
        }
    }
    // last batch was incomplete
    if (n != batch_size) {
        for (int j = 0; j < threads; j++) {
            bam_destroy1(bam_entries[p][j][i + 1]) ;
            bam_entries[p][j][i + 1] = nullptr ;
        }
    }
    lprint({"Loaded", to_string(n), "BAM reads.."});
    return n == batch_size ;
}

