#ifndef FST_HPP
#define FST_HPP

#include <iterator>
#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>
#include <zlib.h>

#include <kseq.h>

#ifdef __GNUC__
#pragma GCC system_header
#endif

KSEQ_INIT(gzFile, gzread)

using namespace std;

struct fastq_entry_t {
  string head;
  string seq;
  string qual;
  uint start;
  uint len;
  // constructor
  fastq_entry_t(const string &h, const string &s, const string &q,
                const uint st = 0, const uint l = 0)
      : head(h), seq(s), qual(q) {
    len = l;
    start = st;
  }
  bool operator==(const fastq_entry_t &o) const { return seq == o.seq; }
};

namespace std {
template <> struct hash<fastq_entry_t> {
  size_t operator()(const fastq_entry_t &k) const {
    return hash<string>()(k.seq);
  }
};
} // namespace std

#endif
