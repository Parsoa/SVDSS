#ifndef BAM_HPP
#define BAM_HPP

#include <mutex>
#include <vector>
#include <string>
#include <utility>
#include <iterator>
#include <unordered_map>

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <htslib/hts_endian.h>

extern uint32_t cigar_len_mask ; 
extern uint32_t cigar_type_mask ;

// CHECKME: I changed this struct but I didn't check if I broke something in realignment
struct CIGAR {
    std::vector<std::pair<uint, char>> ops;
    uint mismatches;
    uint ngaps;
    uint score;

    CIGAR() { 
        mismatches = 0;
        ngaps = -1;
        score = -1;
    }

    CIGAR(std::vector<std::pair<uint, char>> ops_, int score_)
    {
        mismatches = -1;
        score = score_;
        ops = ops_;
        ngaps = 0;
        for (uint i = 0; i < ops_.size(); ++i)
            ngaps += ((ops_[i].second == 'I' || ops_[i].second == 'D') ? 1 : 0);
    }

    void add(int l, char op, int e) {
        mismatches += e;
        if (ops.empty() || ops.back().second != op) {
            ops.push_back(std::make_pair(l, op));
        } else {
            ops.back().first += l;
        }
    }

    void add_front(int l) {
        ops.front().first += l;
        ops.insert(ops.begin(), std::make_pair(1, 'M'));
    }

    void fixclips()
    {
        if (ops.front().second != 'M')
        {
            switch (ops.front().second)
            {
            case 'I':
                ops.front().second = 'S';
                break;
            case 'D':
                ops.erase(ops.begin()); // CHECKME: by doing this we may reduce the length of merged variations
                --ngaps;
                break;
            }
        }
        else if (ops.back().second != 'M')
        {
            switch (ops.back().second)
            {
            case 'I':
                ops.back().second = 'S';
                break;
            case 'D':
                ops.pop_back(); // CHECKME: by doing this we may reduce the length of merged variations
                --ngaps;
                break;
            }
        }
    }

    const std::pair<uint, char> &operator[](std::size_t i) const
    {
        return ops[i];
    }

    uint size() const
    {
        return ops.size();
    }

    std::string to_str() const
    {
        std::string cigar_str;
        for (const auto &op : ops)
            cigar_str += std::to_string(op.first) + op.second;
        return cigar_str;
    }
};

std::string print_cigar_symbol(int type) ;
std::vector<std::pair<uint32_t, uint32_t>> decode_cigar(bam1_t* read) ;
uint8_t* encode_cigar(std::vector<std::pair<uint32_t, uint32_t>> cigar) ;
uint8_t* encode_bam_seq(char* seq) ;
char reverse_complement_base(char base) ;
std::vector<std::pair<int, int>> get_aligned_pairs(bam1_t *alignment) ;

#endif
