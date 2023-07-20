#ifndef BAM_HPP
#define BAM_HPP

#include <ctype.h>
#include <iostream>
#include <iterator>
#include <mutex>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <htslib/hts.h>
#include <htslib/hts_endian.h>
#include <htslib/sam.h>

using namespace std;

extern uint32_t cigar_len_mask;
extern uint32_t cigar_type_mask;

// CHECKME: I changed this struct but I didn't check if I broke something in
// realignment
struct CIGAR {
  vector<pair<uint, char>> ops;
  int mismatches;
  uint ngaps;
  uint start;
  int score;

  CIGAR() {
    mismatches = 0;
    ngaps = -1;
    score = -1;
  }

  CIGAR(vector<pair<uint, char>> ops_, int score_, uint start_ = 0) {
    mismatches = -1;
    score = score_;
    ops = ops_;
    ngaps = 0;
    start = start_;
    for (uint i = 0; i < ops_.size(); ++i) {
      ngaps += ((ops_[i].second == 'I' || ops_[i].second == 'D') ? 1 : 0);
    }
  }

  void parse_cigar(char *cigar) {
    int b = 0;
    int lc = strlen(cigar);
    for (int i = 0; i < lc; i++) {
      if (isdigit(cigar[i])) {
        continue;
      } else {
        char type = cigar[i];
        cigar[i] = '\0';
        int l = stoi(string(cigar + b));
        ops.push_back(make_pair(l, type));
        cigar[i] = type;
        b = i + 1;
      }
    }
  }

  CIGAR(char *cigar, int score_, int start_ = 0) {
    mismatches = -1;
    score = score_;
    start = start_;
    parse_cigar(cigar);
    ngaps = 0;
    for (uint i = 0; i < ops.size(); ++i) {
      ngaps += ((ops[i].second == 'I' || ops[i].second == 'D') ? 1 : 0);
    }
  }

  void add(int l, char op, int e) {
    mismatches += e;
    if (ops.empty() || ops.back().second != op) {
      ops.push_back(make_pair(l, op));
    } else {
      ops.back().first += l;
    }
  }

  void add_front(int l) {
    ops.front().first += l;
    ops.insert(ops.begin(), make_pair(1, 'M'));
  }

  void print() {
    for (uint i = 0; i < ops.size(); ++i) {
      cout << ops[i].first << ops[i].second;
    }
    cout << endl;
  }

  void fixclips() {
    if (ops.front().second != 'M') {
      if (ops.front().second == 'I') {
        ops.front().second = 'S';
      }
      if (ops.front().second == 'D') {
        ops.erase(ops.begin()); // CHECKME: by doing this we may reduce the
                                // length of merged variations
        --ngaps;
      }
    } else if (ops.back().second != 'M') {
      if (ops.back().second == 'I') {
        ops.back().second = 'S';
      }
      if (ops.back().second == 'D') {
        ops.pop_back(); // CHECKME: by doing this we may reduce the length of
                        // merged variations
        --ngaps;
      }
    }
  }

  const pair<uint, char> &operator[](size_t i) const { return ops[i]; }

  uint size() const { return ops.size(); }

  string to_str() const {
    string cigar_str;
    for (const auto &op : ops) {
      cigar_str += to_string(op.first) + op.second;
    }
    return cigar_str;
  }
};

string print_cigar_symbol(int type);
vector<pair<int, int>> decode_cigar(bam1_t *read);
uint8_t *encode_cigar(vector<pair<uint32_t, uint32_t>> cigar);
uint8_t *encode_bam_seq(char *seq);
char reverse_complement_base(char base);
vector<pair<int, int>> get_aligned_pairs(bam1_t *alignment);

#endif
