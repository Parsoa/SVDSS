#include <omp.h>
#include <ctime>
#include <chrono>
#include <string>
#include <vector>
#include <thread>
#include <fstream>
#include <iomanip>
#include <assert.h>
#include <iostream>
#include <pthread.h>
#include <unordered_map>

#include "config.hpp"
#include "ping_pong.hpp"

using namespace std ;

#ifdef DEBUG_MODE
#  define DEBUG(x) x
#  define NEBUG(x)
#else
#  define DEBUG(x)
#  define NEBUG(x) x
#endif

// ============================================================================= \\
// ============================================================================= \\
// ============================================================================= \\

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

//void PingPong::ping_pong_search(rld_t *index, const char* seq, const char* qual, vector<fastq_entry_t>& solutions) {
void PingPong::ping_pong_search(rld_t *index, fastq_entry_t fqe, vector<fastq_entry_t>& solutions) {
    int l = fqe.seq.size() ;
    //int l = strlen(seq) ;
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
        solutions.push_back(get_solution(fqe, begin, sfs_len)) ;
        if (!check_solution(index, fqe.seq.substr(begin, sfs_len))) {
            cerr << "Invalid SFS: " << sfs_len << endl ;
        } ;
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
    cout << "Loaded " << n << " BAM reads.." << endl ;
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
            return true ;
        }
    }
    cout << "Loaded " << n << " FASTQ reads.." << endl ;
    return n != 0 ? true : false ;
}

vector<fastq_entry_t> PingPong::process_batch_fastq(rld_t* index, vector<fastq_entry_t> fastq_entries) {
    vector<fastq_entry_t> solutions ;
    for (const auto fastq_entry : fastq_entries) {
        ping_pong_search(index, fastq_entry, solutions) ;
    }
    return solutions ;
}

vector<fastq_entry_t> PingPong::process_batch_bam(rld_t* index, vector<bam1_t*> bam_entries) {
    vector<fastq_entry_t> solutions ;
    //char* seq = (char*) malloc(10000) ;
    //uint32_t len = 0 ;
    //bam1_t* alignment ; 
    //for (int b = 0; b < bam_entries.size(); b++) {
    //    alignment = bam_entries[b] ;
    //    uint32_t l = alignment->core.l_qseq ; //length of the read
    //    if (l > len) {
    //        if (len > 0) {
    //            free(seq) ;
    //        }
    //        len = l ;
    //        seq = (char*) malloc(l + 1) ;
    //    }
    //    uint8_t *q = bam_get_seq(alignment) ; //quality string
    //    for (int i = 0; i < l; i++){
    //        seq[i] = seq_nt16_str[bam_seqi(q, i)]; //gets nucleotide id and converts them into IUPAC id.
    //    }
    //    seq[l] = '\0' ; // null terminate
    //    ping_pong_search(index, seq, seq, solutions) ;
    //}
    return solutions ;
}

void PingPong::output_batch(void* args) {
    auto c = Configuration::getInstance() ;
    int batch = ((OutputBatchArgs*) args)->batch ;
    string path = c->workdir + "/solution_batch_" + std::to_string(batch - 1) + ".fastq" ;
    cout << "Outputting to " << path << endl ;
    std::ofstream o(path) ;
    for (const auto it : search_solutions[batch - 1]) {
        fastq_entry_t fastq_entry = it.first ;
        o << "@" << fastq_entry.head << ".css" << "_" << fastq_entry.start << ":"
            << fastq_entry.start + fastq_entry.len - 1 << ":" << it.second << endl
            << fastq_entry.seq << endl
            << "+" << endl
            << fastq_entry.qual << endl ;
    }
    search_solutions[batch - 1].clear() ;
    // output read ids
    path = c->workdir + "/read_ids_batch_" + std::to_string(batch - 1) + ".fasta" ;
    cout << "Outputting to " << path << endl ;
    std::ofstream f(path) ;
    for (const auto it : read_ids[batch - 1]) {
        fastq_entry_t fastq_entry = it.first ;
        f << ">" << fastq_entry.seq << endl ;
        for (auto id: it.second) {
            f << id << "$" ;
        }
        f << endl ;
    }
    read_ids[batch - 1].clear() ;
}

int PingPong::search() {
    config = Configuration::getInstance() ;
    // parse arguments
    cout << "Restoring index.." << endl ;
    rld_t *index = rld_restore(config->index.c_str()) ;
    cout << "Done." << endl ;
    int mode = 0 ;
    if (config->fastq != "") {
        cout << "FASTQ input: " << config->fastq << endl ;
        fastq_file = gzopen(config->fastq.c_str(), "r") ;
        fastq_iterator = kseq_init(fastq_file) ;
    } else if (config->bam != "") {
        cout << "BAM input.." << endl ;
        bam_file = hts_open(config->bam.c_str(), "r") ;
        bam_header = sam_hdr_read(bam_file) ; //read header
        mode = 1 ;
    } else {
        cerr << "No input file provided, aborting.." << endl ;
        exit(1) ;
    }
    cout << "Extracting SFS strings.." << endl ;
    cout << "Overlap = " << config->overlap << endl ;
    cout << "Minimum length = " << config->min_string_length << endl ;
    // load first batch
    unordered_map<fastq_entry_t, int> s ;
    search_solutions.push_back(s) ;
    unordered_map<fastq_entry_t, vector<string>> r ;
    read_ids.push_back(r) ;
    vector<vector<vector<fastq_entry_t>>> batches ;
    for(int i = 0; i < 2; i++) {
        bam_entries.push_back(vector<vector<bam1_t*>>(config->threads)) ;
        fastq_entries.push_back(vector<vector<fastq_entry_t>>(config->threads)) ; // current and next output
        batches.push_back(vector<vector<fastq_entry_t>>(config->threads)) ; // previous and current output
    }
    int p = 0 ;
    int batch_size = 10000 ;
    cerr << "Loading first batch" << endl ;
    if (mode == 0) {
        load_batch_fastq(config->threads, batch_size, p) ;
    } else {
        cout << "Allocating BAM buffers.." << endl ;
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < config->threads; j++) {
                for (int k = 0; k <= batch_size / config->threads; k++) {
                    bam_entries[i][j].push_back(bam_init1()) ;
                }
            }
        }
        load_batch_bam(config->threads, batch_size, p) ;
    }
    // main loop
    time_t t ;
    time(&t) ;
    int b = 0 ;
    uint64_t u = 0 ;
    bool loaded_last_batch = false ;
    bool should_load = true ;
    bool should_process = true ;
    while (true) {
        cerr << "Beginning batch " << b + 1 << endl ;
        uint64_t v = u ;
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
                    cerr << "Loaded." << endl ;
                }
            } else if (i == 1) {
                // the pipeline will always merge
                // merge output of previous batch
                if (b >= 1) {
                    int y = 0 ;
                    for (const auto &batch : batches[(p + 1) % 2]) {
                        y += batch.size() ;
                        for (const auto fastq_entry : batch) {
                            if (search_solutions[current_batch].find(fastq_entry) == search_solutions[current_batch].end()) {
                                search_solutions[current_batch][fastq_entry] = 0 ;
                            }
                            search_solutions[current_batch][fastq_entry] += 1 ;
                            read_ids[current_batch][fastq_entry].push_back(fastq_entry.head) ;
                        }
                    }
                    cerr << y << " total sequences." << endl ;
                }
                if (search_solutions[current_batch].size() >= 10000000 || (!should_process && !should_process)) {
                    cerr << "Memory limit reached, dumping output batch " << current_batch << ".." << endl ;
                    current_batch += 1 ;
                    unordered_map<fastq_entry_t, int> s ;
                    search_solutions.push_back(s) ;
                    unordered_map<fastq_entry_t, vector<string>> r ;
                    read_ids.push_back(r) ;
                    OutputBatchArgs* b_args = new OutputBatchArgs() ;
                    b_args->batch = current_batch ;
                    output_batch((void*) b_args) ;
                }
                cerr << "Merged. " << search_solutions[current_batch].size() << " unique sequences." << endl ;
            } else {
                // process current batch
                if (should_process) {
                    batches[p][i - 2] = process_batch_fastq(index, fastq_entries[p][i - 2]) ;
                }
            }
        }
        if (!should_load) {
            cout << "Processed last batch of inputs." << endl ;
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
        cerr << "Processed batch " << std::left << std::setw(10) << b << ". Reads so far " << std::right << std::setw(12) << u << ". Reads per second: " <<  u / (s - t) << ". Time: " << std::setw(8) << std::fixed << s - t << "\n" ;
    }
    //int y = 0 ;
    //for (const auto &batch : batches[(p + 1) % 2]) {
    //    y += batch.size() ;
    //    for (const auto fastq_entry : batch) {
    //        if (search_solutions[current_batch].find(fastq_entry) == search_solutions[current_batch].end()) {
    //            search_solutions[current_batch][fastq_entry] = 0 ;
    //        }
    //        search_solutions[current_batch][fastq_entry] += 1 ;
    //        read_ids[current_batch][fastq_entry].push_back(fastq_entry.head) ;
    //    }
    //}
    //// last batch
    //current_batch += 1 ;
    //OutputBatchArgs* b_args = new OutputBatchArgs() ;
    //b_args->batch = current_batch ;
    //output_batch((void*) b_args) ;
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
        cout << "Appending index: " << c->append << endl ;
        if ((fp = fopen(c->append.c_str(), "rb")) == 0) {
            cerr << "Failed to open file " << c->append << endl ;
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
    cout << "Restoring index.." << endl ;
    rld_t *index = rld_restore(config->index.c_str()) ;
    int l = q.length() ;
    uint8_t *p = (uint8_t*) q.c_str() ;
    seq_char2nt6(l, p);
    bool found_full = backward_search(index, p, l - 1) ;
    bool found_prefix = backward_search(index, p, l - 2) ;
    bool found_suffix = backward_search(index, p + 1, l - 2) ;
    return !found_full & (found_prefix || found_suffix) ;
}
