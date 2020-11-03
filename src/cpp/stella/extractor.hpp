#ifndef EXT_HPP
#define EXT_HPP

#include <mutex>
#include <vector>
#include <string>
#include <iterator>
#include <unordered_map>

#include "counter.hpp"
#include "bed_utils.hpp"

#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

struct fastq_entry_t {
    std::string head ;
    std::string seq ;
    std::string qual ;
    uint start ;
    uint len ;
    // constructor
    fastq_entry_t(const std::string &h, const std::string &s, const std::string &q, const uint st = 0, const uint l = 0) : head(h), seq(s), qual(q) {
        len = l ;
        start = st ;
    }
    bool operator==(const fastq_entry_t &o) const {
        return seq == o.seq ; 
    }
};

namespace std {
    template <> struct hash<fastq_entry_t> {
        std::size_t operator()(const fastq_entry_t& k) const {
            return std::hash<std::string>()(k.seq) ;
        }
    };
}

class Extractor {

public:

    void run() ;

private:

    void load_strings() ;
    void load_reads() ;
    void load_read_ids() ;
    void dump_reads() ;

    std::unordered_map<std::string, Locus> sequences ;
    std::unordered_map<std::string, std::vector<fastq_entry_t>> reads ;
    std::unordered_map<std::string, std::vector<std::string>> read_ids ;

} ;

#endif
