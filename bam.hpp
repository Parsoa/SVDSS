#ifndef BAM_HPP
#define BAM_HPP

#include <mutex>
#include <vector>
#include <string>
#include <iterator>
#include <unordered_map>

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <htslib/hts_endian.h>

extern uint32_t cigar_len_mask ; 
extern uint32_t cigar_type_mask ;

struct CIGAR {
    std::vector<int> ls;
    std::vector<char> ops;
    uint mismatches;

    CIGAR() { 
        mismatches = 0;
    }

    void add(int l, char op, int e) {
        mismatches += e;
        if (ops.empty() || ops.back() != op) {
            ls.push_back(l);
            ops.push_back(op);
        } else {
            ls.back() += l;
        }
    }

    void add_front(int l) {
        ls[0] += l;
        ls.insert(ls.begin(), 1);
        ops.insert(ops.begin(), 'M');
    }

    void fixclips() {
        if (ops.front() != 'M') {
            ops.front() = 'S';
        } 
        if (ops.back() != 'M') {
            ops.back() = 'S';
        }
    }

    std::string to_str() {
        std::string cigarstring;
        for (int i = 0; i < ls.size(); ++i) {
            cigarstring += std::to_string(ls[i]) + ops[i];
        }
        return cigarstring;
    }
};

std::string print_cigar_symbol(int type) ;
std::vector<std::pair<uint32_t, uint32_t>> decode_cigar(bam1_t* read) ;
uint8_t* encode_cigar(std::vector<std::pair<uint32_t, uint32_t>> cigar) ;
uint8_t* encode_bam_seq(char* seq) ;
char reverse_complement_base(char base) ;
std::vector<std::pair<int, int>> get_aligned_pairs(bam1_t *alignment) ;

#endif
