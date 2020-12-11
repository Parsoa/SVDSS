#ifndef EXT_HPP
#define EXT_HPP

#include <mutex>
#include <vector>
#include <string>
#include <iterator>
#include <unordered_map>

#include "fastq.hpp"
#include "counter.hpp"
#include "bed_utils.hpp"

struct SeqReadMap {
    std::string seq ;
    std::string name ;
    std::vector<std::string> read_ids ;
} ;

class Extractor {

public:

    void run() ;

private:

    void load_strings() ;
    void load_reads() ;
    void load_read_ids() ;
    void dump_reads() ;

    std::unordered_map<std::string, SeqReadMap> sequences ;
    //std::unordered_map<std::string, fastq_entry_t> reads ;
    std::unordered_map<std::string, std::vector<std::string>> read_ids ;

} ;

#endif
