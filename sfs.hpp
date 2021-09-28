#ifndef SFS_HPP
#define SFS_HPP

#include <map>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>

static const char RCN[128] = {
    0,   0,   0, 0,   0,   0,   0,   0,   0,   0,   //  0
    0,   0,   0, 0,   0,   0,   0,   0,   0,   0,   // 10
    0,   0,   0, 0,   0,   0,   0,   0,   0,   0,   // 20
    0,   0,   0, 0,   0,   0,   0,   0,   0,   0,   // 30
    0,   0,   0, 0,   0,   0,   0,   0,   0,   0,   // 40
    0,   0,   0, 0,   0,   0,   0,   0,   0,   0,   // 50
    0,   0,   0, 0,   0,   'T', 0,   'G', 0,   0,   // 60
    0,   'C', 0, 0,   0,   0,   0,   0,   'N', 0,   // 70
    0,   0,   0, 0,   'A', 0,   0,   0,   0,   0,   // 80
    0,   0,   0, 0,   0,   0,   0,   'T', 0,   'G', // 90
    0,   0,   0, 'G', 0,   0,   0,   0,   0,   0,   // 100
    'N', 0,   0, 0,   0,   0,   'A', 0,   0,   0,   // 110
    0,   0,   0, 0,   0,   0,   0,   0              // 120
};

struct SFS {
  uint s ;
  uint l ;
  uint c ;
  bool isreversed;

  SFS() {
    s = 0 ;
    l = 0 ;
    c = 0 ;
    isreversed = false;
  }

  SFS(uint s_, uint l_, uint c_, bool isreversed_) {
    s = s_ ;
    l = l_ ;
    c = c_ ;
    isreversed = isreversed_;
  }

  void reverse(uint p) { s = p - s - l; }
};

bool operator<(const SFS &, const SFS &);

struct ExtSFS {
    std::string chrom ;
    std::string qname ;
    int s ;
    int e ;

    ExtSFS(const std::string& _chrom, const std::string& _qname, int _s, int _e) {
        chrom = _chrom ;
        qname = _qname ;
        s = _s ;
        e = _e ;
    }
};

class Consensus {
public:
    std::string seq ;
    std::string chrom ;
    std::string cigar ;
    int s ;
    int e ;

    Consensus(const std::string _seq, const std::string _cigar, const std::string _chrom, int _s, int _e) {
        seq = _seq ;
        cigar = _cigar ;
        chrom = _chrom ;
        s = _s ;
        e = _e ;
    }
};

std::map<std::string, std::vector<SFS>> parse_sfsfile(const std::string &, int);

#endif
