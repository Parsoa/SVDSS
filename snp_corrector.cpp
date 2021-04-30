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

//TODO: why did I have to redefine this?
#define bam_set_seqi(s,i,b) ((s)[(i)>>1] = ((s)[(i)>>1] & (0xf0 >> ((~(i)&1)<<2))) | ((b)<<((~(i)&1)<<2)))

void SnpCorrector::run() {
    // load all variants
    // two pass approach
    // First pass: For each variant check how many of the reads including it actually have one of its alleles
    // Second pass: For each variant with enough coverage, correct all the reads that have it, this can be integrerated into ping-pong directly to save time but for now let's dump a new BAM file
    auto c = Configuration::getInstance() ;
    load_chromosomes(c->reference) ;
    //vcf_variants = load_vcf_file(c->vcf) ;
    //auto _variants = load_vcf_file(c->vcf) ;
    //for (const auto& chrom: _variants) {
    //    vcf_variants[chrom.first] ;
    //    // Keep small INDELs only
    //    std::copy_if (chrom.second.begin(), chrom.second.end(), std::back_inserter(vcf_variants[chrom.first]), [](const vcf_variant_t& v){ return v.svlen <= 5; } );
    //}
    pass(0) ;
}

int binary_search(vector<vcf_variant_t>& a, int begin, int end, int q, int dir) {
    //cout << "Binary search between " << begin << " and " << end << endl ;
    int m = (begin + end) / 2 ;
    if (begin == end || m == begin) {
        return  -1 ;
    }
    if (a[m].pos + a[m].svlen >= q && a[m].pos <= q) {
        // for SNPs this is the as checking a[m].pos == q
        return m ;
    }
    if (a[m].pos < q && a[m].pos + a[m].svlen < q) {
        if (m != a.size() - 1) {
            if (a[m + 1].pos >= q) {
                if (dir == 0) {
                    return m + 1 ;
                } else {
                    return m ;
                }
            }
            return binary_search(a, m, end, q, dir) ;
        } else {
            if (dir == 1) {
                return m ;
            } else {
                return -1 ;
            }
        }
    }
    if (a[m].pos > q) {
        if (m != 0) {
            if (a[m - 1].pos <= q) {
                if (dir == 0) {
                    return m ;
                } else {
                    return m - 1 ;
                }
            }
            return binary_search(a, begin, m, q, dir) ;
        } else {
            if (dir == 0) {
                return m ;
            } else {
                return -1 ;
            }
        }
    }
    assert("Shouldn't be here.") ;
    return -1 ;
}

std::pair<int, int> SnpCorrector::find_variants_in_read(int pos, int len, string chrom) {
    if (vcf_variants.find(chrom) == vcf_variants.end()) {
        return std::make_pair(-1, -1) ;
    }
    // bonary search in variants on the given chromosome
    vector<vcf_variant_t>& chrom_variants = vcf_variants[chrom] ;
    //cout << "Searching for read @" << chrom << " between " << pos << " - " << pos + len << " among " << chrom_variants.size() << " variants." << endl ;
    int begin = binary_search(chrom_variants, 0, chrom_variants.size() - 1, pos, 0) ;
    if (begin == -1) {
        return std::make_pair(-1, -1) ;
    }
    int end = binary_search(chrom_variants, begin + 1, chrom_variants.size() - 1, pos + len - 1, 1) ;
    if (begin != -1 && end != -1) {
        assert(chrom_variants[begin].pos >= pos && chrom_variants[end].pos <= pos + len - 1) ;
    }
    return std::make_pair(begin, end) ;
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

vector<pair<int, int>> calc_cigar_offsets(bam1_t* read) {
    // get CIGAR
    vector<pair<int, int>> cigar_offests ;
    uint32_t* cigar = bam_get_cigar(read) ;
    uint32_t cigar_len_mask = 0xFFFFFFF0 ;
    uint32_t cigar_type_mask = 0xF ;
    int offset = 0 ;
    for (int i = 0; i < read->core.n_cigar; i++) {
        uint32_t type = cigar[i] & cigar_type_mask ;
        uint32_t length = (cigar[i] & cigar_len_mask) >> 4 ;
        cigar_offests.push_back(make_pair(length, type)) ;
        //cout << length << print_cigar_symbol(type) ;
    }
    //cout << endl ;
    return cigar_offests ;
}

//void correct_snps(bam1_t* alignment, std::pair<int, int> limits, char* read_seq, string chrom) {
//    // if read mapped to reverse strand
//    //if ((alignment->core.flag & ((uint16_t) 16)) != 0) {
//        //cout << "Reverse complementing read.." << endl ;
//        //reverse_complement_read(seq) ;
//    //}
//    // finding matching variants
//    //char* neq_seq = (char*) malloc(sizeof(char) * (l + 1)) ;
//    int m = 0 ;
//    int offset = 0 ;
//    int ins_offset = 0 ;
//    int del_offset = 0 ;
//    int group_offset = 0 ; // offset inside the current M group
//    int soft_clip_offset = 0 ;
//    int ref_index = 0 ;
//    int pos = alignment->core.pos + 1 ; // this is 0-based, variant cpoordinates are 1-based
//    // first indel only partially overlap
//    int f = limits.first ;
//    int first_variant_offset = 0 ;
//    auto& first_var = vcf_variants[chrom][f] ;
//    if (first_var.pos < pos) {
//        first_variant_offset = pos - first_var.pos ;
//    }
//    // iterate over all variants in read
//    for (int i = limits.first; i <= limits.second; i++) {
//        if (i != limits.first) {
//            first_variant_offset = 0 ;
//        }
//        auto& var = vcf_variants[chrom][i] ;
//        //cout << "Finding alleles for " << var.pos << endl ;
//        // find position in read corresponding to variant's reference position
//        // This is very tricky. Assume for simplicity that all SNPs are included inside M segments. and so will be sequencing errors.
//        // Now indels may go inside DEL/INS
//        // BAM_CIGAR_TYPE  QUERY  REFERENCE
//        // --------------------------------
//        // BAM_CMATCH      1      1
//        // BAM_CINS        1      0
//        // BAM_CDEL        0      1
//        // BAM_CREF_SKIP   0      1
//        // BAM_CSOFT_CLIP  1      0
//        // BAM_CHARD_CLIP  0      0
//        // BAM_CPAD        0      0
//        // BAM_CEQUAL      1      1
//        // BAM_CDIFF       1      1
//        // BAM_CBACK       0      0
//        // --------------------------------
//        bool fail = false ;
//        while (true) {
//            int diff = var.pos - pos - (offset + del_offset) ;
//            //cout << "Diff: " << diff << ", ref offset: " << offset + del_offset << ", read offset: " << soft_clip_offset + ins_offset + offset << " CIGAR at " << m << "/" << cigar_offests.size() << endl ;
//            if (m == cigar_offests.size()) {
//                //TODO: some late SNPs
//                break ;
//            }
//            if (cigar_offests[m].second == BAM_CMATCH || cigar_offests[m].second == BAM_CEQUAL || cigar_offests[m].second == BAM_CDIFF) {
//                if (diff == 0) {
//                    break ;
//                } else if (diff < 0) {
//                    // overshot, probably because of a large deletion in the read somewhere, can't do this variant anymore
//                    fail = true ;
//                    break ;
//                } else {
//                    if (cigar_offests[m].first - group_offset <= diff) {
//                        // use entire block
//                        strncpy(ref_seq + ref_index, seq + soft_clip_offset + offset + ins_offset, cigar_offests[m].first - group_offset) ;
//                        ref_index += cigar_offests[m].first - group_offset ;
//                        offset += cigar_offests[m].first - group_offset ;
//                    } else {
//                        // use block partially and remain in block for next pass
//                        offset += diff ;
//                        group_offset += diff ;
//                        break ;
//                    }
//                }
//            } else if (cigar_offests[m].second == BAM_CINS) {
//                ins_offset += cigar_offests[m].first ;
//            } else if (cigar_offests[m].second == BAM_CDEL) {
//                del_offset += cigar_offests[m].first ;
//                for (int j = 0; j < cigar_offests[m].first; j++) {
//                    ref_seq[ref_index] = 'N' ;
//                    ref_index++ ;
//                }
//            } else if (cigar_offests[m].second == BAM_CSOFT_CLIP) {
//                if (m == 0) {
//                    soft_clip_offset += cigar_offests[m].first ;
//                }
//            } else if (cigar_offests[m].second == BAM_CREF_SKIP) {
//                // What is even this?
//            } else {//if (cigar_offests[m].second == BAM_CPAD || cigar_offests[m].second == BAM_CHARD_CLIP || cigar_offests[m].second == BAM_CBACK) {
//                // pass
//            }
//            m += 1 ;
//            group_offset = 0 ;
//        }
//        ref_seq[ref_index] = '\0' ;
//        if (fail) {
//            continue ;
//        }
//        int match = -1 ;
//        for (int k = 0; k < 2; k++) {
//            if (var.alleles[k] == "$") {
//                break ;
//            }
//            match = k ;
//            int l = var.alleles[k].length() ;
//            for (int j = 0; j < l; j++) {
//                if (seq[soft_clip_offset + ins_offset + offset + j] != var.alleles[k][j]) {
//                    match = -1 ;
//                    break ;
//                }
//            }
//            if (match != -1) {
//                // rewrite BAM seq
//                // won't change CIGAR and won't change quality scores
//                for (int j = 0; j < l; j++) {
//                    int index = soft_clip_offset + ins_offset + offset + j ;
//                    //cout << "Modified position " << index << " from " << var.alleles[k][j] << " to " << var.ref[j] << endl ;
//                    bam_set_seqi(q, index, seq_nt16_table[(unsigned char) var.ref[j]]) ;
//                }
//                //bam_write1(out_bam_file->fp.bgzf, alignment) ;
//                //should_quit = true ;
//                break ;
//            }
//        }
//        if (match) {
//            //cout << "Matched " << var.pos << "@" << var.chrom << ".." << endl ;
//        }
//        if (!match) {
//            //cout << "Not matched " << var.pos << "@" << var.chrom << ".." << endl ;
//            output[var] += 1 ; // initializes to 0
//        }
//    }
//}

fastq_entry_t SnpCorrector::correct_read(bam1_t* alignment, char* read_seq, string chrom) {
    auto cigar_offests = calc_cigar_offsets(alignment) ;
    int l = 0 ;
    for (auto p: cigar_offests) {
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
    int pos = alignment->core.pos + 1 ; // this is 0-based, variant cpoordinates are 1-based
    while (true) {
        if (m == cigar_offests.size()) {
            break ;
        }
        if (cigar_offests[m].second == BAM_CMATCH || cigar_offests[m].second == BAM_CEQUAL || cigar_offests[m].second == BAM_CDIFF) {
            for (int j = 0; j < cigar_offests[m].first; j++) {
                new_seq[n] = chromosome_seqs[chrom][ref_offset + j] ;
                n++ ;
            }
            ref_offset += cigar_offests[m].first ;
            match_offset += cigar_offests[m].first ;
        } else if (cigar_offests[m].second == BAM_CINS) {
            if (cigar_offests[m].first <= 10) {
                // if a short INDEL then just don't add it to read
            } else {
                // for long INS, this is probably a SV so add it to the read
                for (int j = 0; j < cigar_offests[m].first; j++) {
                    new_seq[n] = read_seq[soft_clip_offset + match_offset + ins_offset + j] ;
                    n++ ;
                }
            }
            ins_offset += cigar_offests[m].first ;
        } else if (cigar_offests[m].second == BAM_CDEL) {
            if (cigar_offests[m].first <= 10) {
                // if a short DEL so let's just fix it
                for (int j = 0; j < cigar_offests[m].first; j++) {
                    new_seq[n] = chromosome_seqs[chrom][ref_offset + j] ;
                    n++ ;
                }
            } else {
                // for long DEL, this is probably a SV so let it be what it was
            }
            del_offset += cigar_offests[m].first ;
            ref_offset += cigar_offests[m].first ;
        } else if (cigar_offests[m].second == BAM_CSOFT_CLIP) {
            for (int j = 0; j < cigar_offests[m].first; j++) {
                new_seq[n] = read_seq[soft_clip_offset + match_offset + ins_offset + j] ;
                n++ ;
            }
            soft_clip_offset += cigar_offests[m].first ;
        } else if (cigar_offests[m].second == BAM_CREF_SKIP) {
            // What is even this?
        } else {//if (cigar_offests[m].second == BAM_CPAD || cigar_offests[m].second == BAM_CHARD_CLIP || cigar_offests[m].second == BAM_CBACK) {
            // pass
        }
        m += 1 ;
    }
    new_seq[n] = '\0' ;
    //cout << read_seq << endl ;
    //cout << new_seq << endl ;
    //char* substr = (char*) malloc(1000) ;
    //strncpy(substr, chromosome_seqs[chrom] + pos - 1, 999) ;
    //substr[999] = '\0' ;
    //cout << substr << endl ; 
    //
    string s(new_seq) ;
    string qname(bam_get_qname(alignment)) ;
    fastq_entry_t f {qname, s, s, pos, n} ;
    free(new_seq) ;
    return f ;
}

vector<fastq_entry_t> SnpCorrector::process_batch_1(vector<bam1_t*> bam_entries) {
    vector<fastq_entry_t> output ;
    char* seq = (char*) malloc(10000) ;
    uint32_t len = 0 ;
    bam1_t* alignment ;
    for (int b = 0; b < bam_entries.size(); b++) {
        alignment = bam_entries[b] ;
        if (alignment == nullptr) {
            break ;
        }

        // skip if unmapped, supplementary or secondary
        if (alignment->core.flag & BAM_FUNMAP || alignment->core.flag & BAM_FSUPPLEMENTARY || alignment->core.flag & BAM_FSECONDARY)
          continue;

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
        if (l < 2) {
            //cerr << "Read too short, ignoring.." << endl ;
            continue ;
        }
        if (alignment->core.tid < 0) {
            continue ;
        }
        string chrom(bam_header->target_name[alignment->core.tid]) ;
        //correct_snps(alignment, limits, seq, chrom) ;
        if (chromosome_seqs.find(chrom) != chromosome_seqs.end()) {
          fastq_entry_t fastq_entry = correct_read(alignment, seq, chrom) ;
          output.push_back(fastq_entry) ;
        }
    }
    free(seq) ;
    return output ;
}

// BAM writing based on https://www.biostars.org/p/181580/
int SnpCorrector::pass(int index) {
    cout << "Running first pass.." << endl ;
    auto config = Configuration::getInstance() ;
    // parse arguments
    bam_file = hts_open(config->bam.c_str(), "r") ;
    bam_header = sam_hdr_read(bam_file) ; //read header
    auto out_path = config->workdir + "/corrected.fastq" ;
    //cout << "Writing correct BAM to " << out_path << endl ;
    //out_bam_file = sam_open(out_path.c_str(), "wb") ;
    //int r = bam_hdr_write(out_bam_file->fp.bgzf, bam_header) ;
    //if (r < 0) {
    //    cerr << "Can't write corrected BAM header, aborting.." << endl ;
    //}
    std::ofstream out_file(out_path) ;
    // confidence scores 
    //vector<vector<unordered_map<vcf_variant_t, int>>> batches ;
    vector<vector<vector<fastq_entry_t>>> batches ;
    for(int i = 0; i < 2; i++) {
        bam_entries.push_back(vector<vector<bam1_t*>>(config->threads)) ;
        //batches.push_back(vector<unordered_map<vcf_variant_t, int>>(config->threads)) ; // previous and current output
        batches.push_back(vector<vector<fastq_entry_t>>(config->threads)) ; // previous and current output
    }
    int p = 0 ;
    int batch_size = 10000 ;
    cerr << "Loading first batch" << endl ;
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
    int b = 0 ;
    bool loaded_last_batch = false ;
    bool should_load = true ;
    bool should_process = true ;
    uint64_t u = 0 ;
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
                    int ret ;
                    for (int j = 0; j < config->threads; j++) {
                        for (auto& fastq_entry: batches[(p + 1) % 2][j]) {
                            out_file << "@" << fastq_entry.head << endl
                                << fastq_entry.seq << endl
                                << "+" << endl
                                << fastq_entry.seq << endl ;
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
            } else if (i == 1) {
                // merge output of previous batch
                // nothing to be merged
                //if (b >= 1) {
                //    int y = 0 ;
                //    for (const auto &batch : batches[(p + 1) % 2]) {
                //        y += batch.size() ;
                //       for (auto it = batch.begin(); it != batch.end() ;it++) {
                //            confidence_scores[it->first] += it->second ;
                //        }
                //    }
                //    cerr << y << " total sequences." << endl ;
                //}
                //cerr << "Merged. " << confidence_scores.size() << " unique sequences." << endl ;
            } else {
                // process current batch
                if (should_process) {
                    batches[p][i - 2] = process_batch_1(bam_entries[p][i - 2]) ;
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
    cout << "Done." << endl ;
    sam_close(bam_file) ;
    //sam_close(out_bam_file) ;
    return 0 ;
}

bool SnpCorrector::load_batch_bam(int threads, int batch_size, int p) {
    int n = 0 ;
    int i = 0 ;
    while (sam_read1(bam_file, bam_header, bam_entries[p][n % threads][i]) >= 0) {
        n += 1 ;
        if (n % threads == threads - 1) {
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

void write_bam_entry(bam1_t* b) {

}

void SnpCorrector::stats() {
    int n = 0;
    for (auto it = vcf_variants.begin(); it != vcf_variants.end(); it++) {
        n += it->second.size() ;
    }
    cout << n << " total variants." << ends ;
    cout << confidence_scores.size() << " variants with missing reads.." << endl ;
}
