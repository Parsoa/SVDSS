#include <omp.h>
#include <ctime>
#include <chrono>
#include <thread>
#include <fstream>
#include <iomanip>
#include <assert.h>
#include <iostream>
#include <pthread.h>

#include "config.hpp"
#include "chromosomes.hpp"
#include "snp_corrector.hpp"

using namespace std ;

#define BAM_CIGAR_MASK  0xf
#define BAM_CIGAR_TYPE  0x3C1A7

void SnpCorrector::run() {
    auto c = Configuration::getInstance() ;
    load_chromosomes(c->reference) ;
    correct_reads() ;
}

char reverse_complement_base(char base) {
    if (base == 'C' || base == 'c') {
        return 'G' ;
    }
    if (base == 'A' || base == 'a') {
        return 'T' ;
    }
    if (base == 'G' || base == 'g') {
        return 'C' ;
    }
    if (base == 'T' || base == 't') {
        return 'A' ;
    }
    else {
        return 'N' ;
    }
}

void reverse_complement_read(char* seq) {
    int l = strlen(seq) ;
    int i = 0 ;
    while (i < l / 2) {
        auto t = reverse_complement_base(seq[l - i]) ;
        seq[l - 1 - i] = reverse_complement_base(seq[i]) ;
        seq[i] = t ;
        i += 1 ;
    }
}

string print_cigar_symbol(int type) {
    if (type == BAM_CMATCH) {
        return "M" ;
    }
    if (type == BAM_CINS) {
        return "I" ;
    }
    if (type == BAM_CDEL) {
        return "D" ;
    }
    if (type == BAM_CSOFT_CLIP) {
        return "S" ;
    }
    if (type == BAM_CHARD_CLIP) {
        return "H" ;
    }
    return "X" ;
}

uint32_t cigar_len_mask = 0xFFFFFFF0 ;
uint32_t cigar_type_mask = 0xF ;

vector<pair<uint32_t, uint32_t>> decode_cigar(bam1_t* read) {
    // get CIGAR
    vector<pair<uint32_t, uint32_t>> cigar_offsets ;
    uint32_t* cigar = bam_get_cigar(read) ;
    int offset = 0 ;
    for (int i = 0; i < read->core.n_cigar; i++) {
        uint32_t type = cigar[i] & cigar_type_mask ;
        uint32_t length = cigar[i] >> 4 ;
        cigar_offsets.push_back(make_pair(length, type)) ;
    }
    return cigar_offsets ;
}

uint8_t* encode_cigar(vector<pair<uint32_t, uint32_t>> cigar) {
    uint32_t* cigar_bytes = (uint32_t*) malloc(sizeof(uint32_t) * cigar.size()) ;
    for (int i = 0; i < cigar.size(); i++) {
        cigar_bytes[i] = (cigar[i].first << 4) | (cigar[i].second & cigar_type_mask) ; 
    }
    return (uint8_t*) cigar_bytes ;
}

uint8_t* encode_bam_seq(char* seq) {
    int n = (strlen(seq) + 1) >> 1 ;
    int l_seq = strlen(seq) ;
    //cout << "l_seq: " << l_seq << " " << n << endl ;
    uint8_t* seq_bytes = (uint8_t*) malloc(sizeof(uint8_t) * n) ;
    int i = 0 ;
    n = 0 ;
    for (i = 0; i + 1 < l_seq; i += 2) {
        seq_bytes[n] = (seq_nt16_table[(unsigned char)seq[i]] << 4) | seq_nt16_table[(unsigned char)seq[i + 1]];
        n += 1 ;
    }
    for (; i < l_seq; i++) {
        seq_bytes[n] = seq_nt16_table[(unsigned char)seq[i]] << 4;
        n += 1 ;
    }
    //cout << "i " << i << " " << n << endl ;
    return seq_bytes ;
}

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
//    uint16_t bin;       // TODO 
//    uint8_t qual;       // won't change
//    uint8_t l_extranul; // won't change 
//    uint16_t flag;      // won't change
//    uint16_t l_qname;   // won't change
//    uint32_t n_cigar;   // will change
//    int32_t l_qseq;     // will change
//    int32_t mtid;       // TODO won't change
//    hts_pos_t mpos;     // TODO won't change
//    hts_pos_t isize;    // may or may ont change?
//} bam1_core_t;

void rebuild_bam_entry(bam1_t* alignment, char* seq, uint8_t* qual, vector<pair<uint32_t, uint32_t>> cigar) {
    int data_size = 0 ;
    // update core
    alignment->core.n_cigar = cigar.size() ;
    alignment->core.l_qseq = strlen(seq) ;
    int l = strlen(seq) ;
    // rebuild data
    alignment->l_data = alignment->core.l_qname + (4 * alignment->core.n_cigar) + ((l + 1) >> 1) + l ; // bam_get_l_aux(alignment) ;
    //uint8_t* new_data = (uint8_t*) malloc(sizeof(uint8_t) * alignment->l_data) ;
    //sam_realloc_bam_data(alignment, alignment->l_data) ;
    //
    // copy qname
    int offset = alignment->core.l_qname ;
    //memcpy(new_data, alignment->data, offset) ;
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
    // don't copy aux
    // memcpy(new_data + offset, bam_get_aux(alignment), bam_get_l_aux(alignment)) ;
    //
    //free(alignment->data) ;
    //alignment->data = new_data ;
}

fastq_entry_t SnpCorrector::correct_read(bam1_t* alignment, char* read_seq, string chrom) {
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
    while (true) {
        if (m == cigar_offsets.size()) {
            break ;
        }
        if (cigar_offsets[m].second == BAM_CMATCH || cigar_offsets[m].second == BAM_CEQUAL || cigar_offsets[m].second == BAM_CDIFF) {
            for (int j = 0; j < cigar_offsets[m].first; j++) {
                new_seq[n] = chromosome_seqs[chrom][ref_offset + j] ;
                new_qual[n] = qual[soft_clip_offset + match_offset + ins_offset + j] ;
                n++ ;
            }
            ref_offset += cigar_offsets[m].first ;
            match_offset += cigar_offsets[m].first ;
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
            // What is even this?
        } else {//if (cigar_offsets[m].second == BAM_CPAD || cigar_offsets[m].second == BAM_CHARD_CLIP || cigar_offsets[m].second == BAM_CBACK) {
            // pass
        }
        m += 1 ;
    }
    new_seq[n] = '\0' ;
    new_qual[n] = '\0' ;
    //cout << "old seq length: " << strlen(read_seq) << " new seq len: " << n << endl ;
    //cout << "old CIGAR length: " << cigar_offsets.size() << " new CIGAR len: " << new_cigar.size() << endl ;
    rebuild_bam_entry(alignment, new_seq, new_qual, new_cigar) ;
    string s(new_seq) ;
    string qname(bam_get_qname(alignment)) ;
    fastq_entry_t f {qname, s, s, pos, n} ;
    free(new_seq) ;
    return f ;
}

vector<fastq_entry_t> SnpCorrector::process_batch(vector<bam1_t*> bam_entries) {
    vector<fastq_entry_t> output ;
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
        //correct_snps(alignment, limits, seq, chrom) ;
        //cout << bam_get_qname(alignment) << " " << bam_header->target_name[alignment->core.tid] << " " << alignment->core.mpos << endl ;
        fastq_entry_t fastq_entry = correct_read(alignment, seq, chrom) ;
        output.push_back(fastq_entry) ;
    }
    free(seq) ;
    return output ;
}

// BAM writing based on https://www.biostars.org/p/181580/
int SnpCorrector::correct_reads() {
    cout << "Running first pass.." << endl ;
    auto config = Configuration::getInstance() ;
    // parse arguments
    bam_file = hts_open(config->bam.c_str(), "r") ;
    bam_header = sam_hdr_read(bam_file) ; //read header
    //auto out_path = config->workdir + "/reconstructed.fastq" ;
    auto out_bam_path = config->workdir + "/reconstructed.bam" ;
    //cout << "Writing correct BAM to " << out_path << endl ;
    out_bam_file = sam_open(out_bam_path.c_str(), "wb") ;
    int r = bam_hdr_write(out_bam_file->fp.bgzf, bam_header) ;
    if (r < 0) {
        cerr << "Can't write reconstructed BAM header, aborting.." << endl ;
    }
    //std::ofstream out_file(out_path) ;
    // confidence scores 
    vector<vector<vector<fastq_entry_t>>> batches ;
    for(int i = 0; i < 2; i++) {
        bam_entries.push_back(vector<vector<bam1_t*>>(config->threads)) ;
        batches.push_back(vector<vector<fastq_entry_t>>(config->threads)) ; // previous and current output
    }
    int p = 0 ;
    int b = 0 ;
    int batch_size = 10000 ;
    cerr << "Loading first batch.." << endl ;
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
        cerr << "Beginning batch " << b + 1 << endl ;
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
                    // write FASTQ
                    //for (int j = 0; j < config->threads; j++) {
                    //    for (auto& fastq_entry: batches[(p + 1) % 2][j]) {
                    //        out_file << "@" << fastq_entry.head << endl
                    //            << fastq_entry.seq << endl
                    //            << "+" << endl
                    //            << fastq_entry.seq << endl ;
                    //    }
                    //}
                    // write BAM
                    int ret = 0 ;
                    for (int k = 0; k < batch_size / config->threads; k++) {
                        for (int j = 0; j < config->threads; j++) {
                            if (bam_entries[(p + 1) % 2][j][k] != nullptr) {
                                auto alignment = bam_entries[(p + 1) % 2][j][k] ;
                                //if (alignment->core.flag & BAM_FUNMAP || alignment->core.flag & BAM_FSUPPLEMENTARY || alignment->core.flag & BAM_FSECONDARY) {
                                //    continue ;
                                //}
                                ret = bam_write1(out_bam_file->fp.bgzf, bam_entries[(p + 1) % 2][j][k]);
                                num_reads++ ;
                                if (ret < 0) {
                                    cerr << "Can't write corrected BAM record, aborting.." << endl ;
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
                        cout << "Last input batch loaded." << endl ;
                    } else {
                        cerr << "Loaded." << endl ;
                    }
                }
                //loaded_last_batch = true ;
            } else if (i == 1) {
                // merge output of previous batch
            } else {
                // process current batch
                if (should_process) {
                    batches[p][i - 2] = process_batch(bam_entries[p][i - 2]) ;
                }
            }
        }
        if (should_terminate) {
            cerr << "Something went wrong, aborting.." << endl ;
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
    cout << "Done." << endl ;
    sam_close(bam_file) ;
    sam_close(out_bam_file) ;
    cout << "Wrote " << num_reads << " reads." << endl ;
    return 0 ;
}

bool SnpCorrector::load_batch_bam(int threads, int batch_size, int p) {
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
    cout << "Loaded " << n << " BAM reads.." << endl ;
    return n == batch_size ;
}

