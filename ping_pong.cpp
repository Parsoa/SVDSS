#include "ping_pong.hpp"

std::string interval2str(rldintv_t sai) {
    return "[" + std::to_string(sai.x[0]) + "," + std::to_string(sai.x[1]) + "," + std::to_string(sai.x[2]) + "]" ;
}

static inline int kputsn(const char *p, int l, kstring_t *s) {
    if (s->l + l + 1 >= s->m) {
        char *tmp;
        s->m = s->l + l + 2;
        kroundup32(s->m);
        if ((tmp = (char*)realloc(s->s, s->m))) s->s = tmp;
        else return EOF;
    }
    memcpy(s->s + s->l, p, l);
    s->l += l;
    s->s[s->l] = 0;
    return l;
}

void seq_char2nt6(int l, unsigned char *s) {
    int i ;
    for (i = 0; i < l; ++i) {
        s[i] = s[i] < 128? seq_nt6_table[s[i]] : 5 ;
    }
}

bool PingPong::check_solution(rld_t* index, std::string S) {
    int l = S.length() ;
    uint8_t *P = (uint8_t*) S.c_str() ;
    seq_char2nt6(l, P);
    bool found_full = backward_search(index, P, l - 1) ;
    bool found_prefix = backward_search(index, P, l - 2) ;
    bool found_suffix = backward_search(index, P + 1, l - 2) ;
    return !found_full & (found_prefix || found_suffix) ;
}

fastq_entry_t PingPong::get_solution(fastq_entry_t fqe, int s, int l) {
    std::string S (fqe.seq, s, l) ;
    std::string Q (fqe.qual, s, l) ;
    return fastq_entry_t(fqe.head, S, Q, s, l) ;
}

bool PingPong::backward_search(rld_t *index, const uint8_t *P, int p2) {
    rldintv_t sai ; // rldintv_t is the struct used to store a SA interval.
    fm6_set_intv(index, P[p2], sai) ;
    while(sai.x[2] != 0 && p2 > 0) {
        --p2;
        rldintv_t osai[6] ;
        rld_extend(index, &sai, osai, 1) ; //1: backward, 0: forward
        sai = osai[P[p2]] ;
    }
    return sai.x[2] != 0 ;
}

void PingPong::ping_pong_search(rld_t *index, fastq_entry_t fqe, std::vector<sfs_solution_t>& solutions) {
    int l = fqe.seq.size() ;
    if (l <= 10) {
        return ;
    }
    DEBUG(cerr << "Read Length: " << l << endl ;)
    char *seq = new char[l + 1] ; // current sequence
    strcpy(seq, fqe.seq.c_str()) ; // seq
    uint8_t *P = (uint8_t*) seq ; 
    seq_char2nt6(l, P) ; // convert to integers
    rldintv_t sai ;

    int begin = l - 1 ;
    while (begin >= 0) {
        // Backward search. Find a mismatching sequence. Stop at first mismatch.
        int bmatches = 0 ;
        fm6_set_intv(index, P[begin], sai) ;
        DEBUG(cerr << "BS from " << int2char[P[begin]] << " (" << begin << "): " << interval2str(sai) << endl ;)
        bmatches = 0 ;
        while (sai.x[2] != 0 && begin > 0) {
            begin-- ;
            bmatches++ ;
            rldintv_t osai[6] ; // output SA intervals (one for each symbol between 0 and 5)
            rld_extend(index, &sai, osai, 1) ;
            sai = osai[P[begin]] ;
            DEBUG(cerr << "- BE with " << int2char[P[begin]] << " (" << begin << "): " << interval2str(sai) << endl ;)
        }
        //last sequence was a match
        if (begin == 0 && sai.x[2] != 0) {
            break ;
        }
        DEBUG(cerr << "Mismatch " << int2char[P[begin]] << " (" <<  begin << "). bmatches: " << to_string(bmatches) << endl ;)
        // Forward search: 
        int end = begin ;
        int fmatches = 0 ;
        fm6_set_intv(index, P[end], sai) ;
        DEBUG(cerr << "FS from " << int2char[P[end]] << " (" << end << "): " << interval2str(sai) << endl ;)
        while(sai.x[2] != 0) {
            end++ ;
            fmatches++ ;
            rldintv_t osai[6] ;
            rld_extend(index, &sai, osai, 0) ;
            sai = osai[P[end] >= 1 && P[end] <= 4 ? 5 - P[end] : P[end]];
            DEBUG(cerr << "- FE with " << int2char[P[end]] << " (" <<  end << "): " << interval2str(sai) << endl ;)
        }
        DEBUG(cerr << "Mismatch " << int2char[P[end]] << " (" << end << "). fmatches: " << fmatches << endl ;)
        // add solution
        DEBUG(cerr << "Adding [" << begin << ", " << end << "]." << endl ;)
        int acc_len = end - begin + 1 ;
        int sfs_len = end - begin + 1 ;
        //if (config->min_string_length > 0) {
        //    sfs_len = acc_len > config->min_string_length ? acc_len : config->min_string_length ;
        //    if (begin + sfs_len >= l - 1) {
        //        sfs_len = acc_len ;
        //    }
        //    assert(sfs_len == config->min_string_length || sfs_len == acc_len) ;
        //    DEBUG(cerr << "Adjusted length to " << sfs_len << "." << endl ;)
        //}
        DEBUG(cerr << "Adjusted length from " << acc_len << " to " << sfs_len << "." << endl ;)

        // CHECKMERGE
        solutions.push_back(sfs_solution_t{begin, sfs_len, fqe.seq.substr(begin, sfs_len)});
        // if (!check_solution(index, fqe.seq.substr(begin, sfs_len))) {
        //   cerr << "Invalid SFS: " << sfs_len << endl ;
        // } ;
        DEBUG(std::this_thread::sleep_for(std::chrono::seconds(1)) ;)

        // prepare for next round
        if (begin == 0) {
            break ;
        }
        if (config->overlap == 0) { // Relaxed
            begin -= 1 ;
        } else {
            if (config->overlap > 0) {
                int overlap = config->overlap >= acc_len ? acc_len - 1 : config->overlap ;
                begin = begin + overlap ;
                assert(begin <= end && overlap >= 0) ;
            } else {
                begin = end + config->overlap ; // overlap < 0
            }
        }
    }
    DEBUG(std::this_thread::sleep_for(std::chrono::seconds(2)) ;)
    delete[] seq ;
}

bool PingPong::load_batch_bam(int threads, int batch_size, int p) {
    for (int i = 0; i < threads; i++) {
        fastq_entries[p][i].clear() ;
    }
    int i = 0 ;
    int n = 0 ;
    while (sam_read1(bam_file, bam_header, bam_entries[p][n % threads][i]) >= 0) {
        auto alignment = bam_entries[p][n % threads][i] ;
        if (alignment == nullptr) {
            break ;
        }
        // these reads have not been reconstructed
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
        uint32_t l = alignment->core.l_qseq ; //length of the read
        char* seq = (char*) malloc(l + 1) ;
        uint8_t *q = bam_get_seq(alignment) ; //quality string
        for (int i = 0; i < l; i++){
            seq[i] = seq_nt16_str[bam_seqi(q, i)]; //gets nucleotide id and converts them into IUPAC id.
        }
        seq[l] = '\0' ; // null terminate
        fastq_entries[p][n % threads].push_back(fastq_entry_t(bam_get_qname(alignment), seq, seq)) ;
        n += 1 ;
        if (n % threads == threads - 1) {
            i += 1 ;
        }
        if (n == batch_size) {
            return true ;
        }
    }
    lprint({"Loaded", to_string(n), "BAM reads.."});
    return n != 0 ? true : false ;
}

bool PingPong::load_batch_fastq(int threads, int batch_size, int p) {
    for (int i = 0; i < threads; i++) {
        fastq_entries[p][i].clear() ;
    }
    int l = 0 ;
    int n = 0 ;
    while ((l = kseq_read(fastq_iterator)) >= 0) {
        fastq_entries[p][n % threads].push_back(fastq_entry_t(fastq_iterator->name.s, fastq_iterator->seq.s, fastq_iterator->qual.s)) ;
        n += 1 ;
        if (n == batch_size) {
          lprint({"Loaded", to_string(n), "FASTQ reads.."});
            return true ;
        }
    }
    lprint({"Loaded", to_string(n), "FASTQ reads.."});
    return n != 0 ? true : false ;
}

batch_type_t PingPong::process_batch(rld_t* index, vector<fastq_entry_t> fastq_entries) {
    batch_type_t solutions ;
    // store read id once for all strings to save space, is it worth it?
    for (const auto fastq_entry: fastq_entries) {
        ping_pong_search(index, fastq_entry, solutions[fastq_entry.head]) ;
    }
    return solutions ;
}

void PingPong::output_batch(int b) {
    auto c = Configuration::getInstance();
    string path = c->workdir + "/solution_batch_" + std::to_string(current_batch) + ".sfs";
    lprint({"Outputting to", path});
    std::ofstream o(path);
    for (int i = last_dumped_batch; i < b; i++) { // for each of the unmerged batches
        for (auto &batch: batches[i]) { // for each thread in batch
            for (auto &read: batch) { // for each read in thread
                int j = 0 ;
                for (auto &sfs: read.second) { // for each sfs in read
                    // optimize file output size by not outputing read name for every SFS
                    o << (j == 0 ? read.first : "*") << "\t" << sfs.seq << "\t" << sfs.begin << "\t" << sfs.len << "\t" << 1 << endl ;
                    j += 1 ;
                }
            }
            batch.clear() ;
        }
        batches[i].clear() ;
    }
    last_dumped_batch = b ;
}

int PingPong::search() {
    config = Configuration::getInstance() ;
    // parse arguments
    lprint({"Restoring index.."});
    rld_t *index = rld_restore(config->index.c_str()) ;
    lprint({"Done."});
    int mode = 0 ;
    if (config->fastq != "") {
        lprint({"FASTQ input:", config->fastq});
        fastq_file = gzopen(config->fastq.c_str(), "r") ;
        fastq_iterator = kseq_init(fastq_file) ;
    } else if (config->bam != "") {
        lprint({"BAM input.."});
        bam_file = hts_open(config->bam.c_str(), "r") ;
        bam_header = sam_hdr_read(bam_file) ; //read header
        mode = 1 ;
    } else {
        lprint({"No input file provided, aborting.."}, 2);
        exit(1) ;
    }
    lprint({"Overlap =", to_string(config->overlap)});
    lprint({"Minimum length =", to_string(config->min_string_length)});
    // load first batch
    for(int i = 0; i < 2; i++) {
        bam_entries.push_back(vector<vector<bam1_t*>>(config->threads)) ;
        fastq_entries.push_back(vector<vector<fastq_entry_t>>(config->threads)) ; // current and next output
    }
    int p = 0 ;
    int batch_size = 10000 ;
    lprint({"Extracting SFS strings on", to_string(config->threads), "threads.."});
    lprint({"Loading first batch.."});
    if (mode == 0) {
        load_batch_fastq(config->threads, batch_size, p) ;
    } else {
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < config->threads; j++) {
                for (int k = 0; k <= batch_size / config->threads; k++) {
                    bam_entries[i][j].push_back(bam_init1()) ;
                }
            }
        }
        load_batch_bam(config->threads, batch_size, p) ;
    }
    batches.push_back(vector<batch_type_t>(config->threads)) ; // previous and current output
    // main loop
    time_t t ;
    time(&t) ;
    int b = 0 ;
    uint64_t u = 0 ;

    bool should_load = true ;
    bool should_process = true ;
    bool loaded_last_batch = false ;
    bool should_update_current_batch = false ;

    int total_sfs = 0 ;
    int total_sfs_output_batch = 0 ;

    while (true) {
        lprint({"Beginning batch", to_string(b + 1)});
        for (int i = 0 ; i < config->threads ; i++) {
            u += fastq_entries[p][i].size() ;
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
                // load next batch of entries
                if (should_load) {
                    if (mode == 0) {
                        loaded_last_batch = !load_batch_fastq(config->threads, batch_size, (p + 1) % 2) ;
                    } else {
                        loaded_last_batch = !load_batch_bam(config->threads, batch_size, (p + 1) % 2) ;
                    }
                    lprint({"Loaded."});
                    batches.push_back(vector<batch_type_t>(config->threads)) ; // previous and current output
                }
            } else if (i == 1) {
               if (b >= 1) {
                    // just count how many strings we have
                    int c = 0 ;
                    for (const auto &batch: batches[b - 1]) {
                        for (auto it = batch.begin(); it != batch.end(); it++) {
                            c += it->second.size() ;
                        } 
                    }
                    total_sfs += c ;
                    total_sfs_output_batch += c ;
                    cerr << "Merged " << c << " new sequences. " << total_sfs << " total sequences." << endl ;
                    // reached memory limit or last pipeline run
                    if (total_sfs_output_batch >= 10000000 || !should_process) {
                        output_batch(b) ;
                        total_sfs_output_batch = 0 ;
                        current_batch += 1 ;
                    }
                }
            } else {
                // process current batch
                if (should_process) {
                    batches[b][i - 2] = process_batch(index, fastq_entries[p][i - 2]) ;
                }
            }
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
    }

    // cleanup
    kseq_destroy(fastq_iterator) ;
    gzclose(fastq_file) ;
    num_output_batches = current_batch ;
    return u ;
}

// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\

/** Code adapted from ropebwt2 (main_ropebwt2 in main.c) **/
int PingPong::index() {
    auto c = Configuration::getInstance() ;
    // hardcoded parameters
    uint64_t m = (uint64_t)(.97 * 10 * 1024 * 1024 * 1024) + 1 ; // batch size for multi-string indexing
    int block_len = ROPE_DEF_BLOCK_LEN, max_nodes = ROPE_DEF_MAX_NODES, so = MR_SO_RCLO ;
    int thr_min = 100 ; // switch to single thread when < 100 strings remain in a batch

    // the index
    mrope_t *mr = 0 ;

    bool binary_output = c->binary ;
    if (c->append != "") {
        FILE *fp ;
        lprint({"Appending index:", c->append});
        if ((fp = fopen(c->append.c_str(), "rb")) == 0) {
            lprint({"Failed to open file", c->append}, 2);
            return 1 ;
        }
        mr = mr_restore(fp) ;
        fclose(fp) ;
    }

    // Initialize mr if not restored
    if (mr == 0) mr = mr_init(max_nodes, block_len, so) ;
    mr_thr_min(mr, thr_min) ;

    // Parsing the input sample
    gzFile fp = gzopen(c->fastq.c_str(), "rb") ;
    kseq_t *ks = kseq_init(fp) ;
    kstring_t buf = { 0, 0, 0 } ; // buffer, will contain the concatenation
    int l ;
    uint8_t *s ;
    int i ;
    while ((l = kseq_read(ks)) >= 0) {
        s = (uint8_t*)ks->seq.s ;

        // change encoding
        for (i = 0; i < l; ++i) {
            s[i] = s[i] < 128? seq_nt6_table[s[i]] : 5 ;
        }

        // Reverse the sequence
        for (i = 0; i < l>>1; ++i) {
            int tmp = s[l-1-i] ;
            s[l-1-i] = s[i] ;
            s[i] = tmp ;
        }

        // Add forward to buffer
        kputsn((char*)ks->seq.s, ks->seq.l + 1, &buf) ;

        // Add reverse to buffer
        for (i = 0; i < l>>1; ++i) {
            int tmp = s[l - 1 - i] ;
            tmp = (tmp >= 1 && tmp <= 4)? 5 - tmp : tmp ;
            s[l - 1 - i] = (s[i] >= 1 && s[i] <= 4)? 5 - s[i] : s[i] ;
            s[i] = tmp ;
        }
        if (l&1) s[i] = (s[i] >= 1 && s[i] <= 4)? 5 - s[i] : s[i] ;
        kputsn((char*)ks->seq.s, ks->seq.l + 1, &buf) ;

        if(buf.l >= m) {
            mr_insert_multi(mr, buf.l, (uint8_t*)buf.s, 1) ;
            buf.l = 0 ;
        }
    }

    if (buf.l) { // last batch
        mr_insert_multi(mr, buf.l, (uint8_t*)buf.s, 1) ;
    }

    free(buf.s) ;
    kseq_destroy(ks) ;
    gzclose(fp) ;

    // dump index to stdout
    if (binary_output) {
        // binary FMR format
        mr_dump(mr, fopen(c->index.c_str(), "wb")) ;
    } else {
        // FMD format
        mritr_t itr ;
        const uint8_t *block ;
        rld_t *e = 0 ;
        rlditr_t di ;
        e = rld_init(6, 3) ;
        rld_itr_init(e, &di, 0) ;
        mr_itr_first(mr, &itr, 1) ;
        while ((block = mr_itr_next_block(&itr)) != 0) {
            const uint8_t *q = block + 2, *end = block + 2 + *rle_nptr(block) ;
            while (q < end) {
                int c = 0 ;
                int64_t l ;
                rle_dec1(q, c, l) ;
                rld_enc(e, &di, l, c) ;
            }
        }
        rld_enc_finish(e, &di) ;
        rld_dump(e, c->index.c_str()) ;
    }

    mr_destroy(mr) ;

    return 0 ;
}

bool PingPong::query(string q) {
    config = Configuration::getInstance() ;
    // parse arguments
    lprint({"Restoring index.."});
    rld_t *index = rld_restore(config->index.c_str()) ;
    int l = q.length() ;
    uint8_t *p = (uint8_t*) q.c_str() ;
    seq_char2nt6(l, p);
    bool found_full = backward_search(index, p, l - 1) ;
    bool found_prefix = backward_search(index, p, l - 2) ;
    bool found_suffix = backward_search(index, p + 1, l - 2) ;
    return !found_full & (found_prefix || found_suffix) ;
}
