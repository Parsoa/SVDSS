#include "bam.hpp"

using namespace std ;

uint32_t cigar_len_mask = 0xFFFFFFF0 ;
uint32_t cigar_type_mask = 0xF ;

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
    return seq_bytes ;
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

vector<pair<int, int>> get_aligned_pairs(bam1_t *alignment) {
    vector<pair<int, int>> result ;
    uint ref_pos = alignment->core.pos ;
    uint read_pos = 0 ;
    auto cigar_offsets = decode_cigar(alignment) ;
    int m = 0 ;
    while (true) {
        if (m == cigar_offsets.size()) {
            break ;
        }
        if (cigar_offsets[m].second == BAM_CMATCH or cigar_offsets[m].second == BAM_CEQUAL or cigar_offsets[m].second == BAM_CDIFF) {
            for (uint i = ref_pos; i < ref_pos + cigar_offsets[m].first; ++i) {
                result.push_back(make_pair(read_pos, i));
                read_pos++;
            }
            ref_pos += cigar_offsets[m].first;
        } else if (cigar_offsets[m].second == BAM_CINS or cigar_offsets[m].second == BAM_CSOFT_CLIP) {
            for (uint i = 0; i < cigar_offsets[m].first; ++i) {
                result.push_back(make_pair(read_pos, -1));
                read_pos++;
            }
        } else if (cigar_offsets[m].second == BAM_CDEL) {
            for (uint i = ref_pos; i < ref_pos + cigar_offsets[m].first; ++i) {
                result.push_back(make_pair(-1, i));
            }
            ref_pos += cigar_offsets[m].first;
        } else if (cigar_offsets[m].second == BAM_CHARD_CLIP) {
            // advances neither
        } else if (cigar_offsets[m].second == BAM_CREF_SKIP) {
            for (uint i = ref_pos; i < ref_pos + cigar_offsets[m].first; ++i) {
                result.push_back(make_pair(-1, i));
            }
            ref_pos += cigar_offsets[m].first;
        } else { //if (cigar_offsets[m].second == BAM_CPAD) {
            //TODO
        }
        m++ ;
    }
    return result;
}

