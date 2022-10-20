#ifndef FST_HPP
#define FST_HPP

#include <iterator>
#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>

#ifdef __GNUC__
#pragma GCC system_header
#endif

#include "kseq.h"
#include <zlib.h>

KSEQ_INIT(gzFile, gzread)

struct fastq_entry_t {
  std::string head;
  std::string seq;
  std::string qual;
  uint start;
  uint len;
  // constructor
  fastq_entry_t(const std::string &h, const std::string &s,
                const std::string &q, const uint st = 0, const uint l = 0)
      : head(h), seq(s), qual(q) {
    len = l;
    start = st;
  }
  bool operator==(const fastq_entry_t &o) const { return seq == o.seq; }
};

namespace std {
template <> struct hash<fastq_entry_t> {
  std::size_t operator()(const fastq_entry_t &k) const {
    return std::hash<std::string>()(k.seq);
  }
};
} // namespace std

#endif
