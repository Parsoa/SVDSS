#ifndef CLUSTER_HPP
#define CLUSTER_HPP

#include <string>
#include <vector>
#include <algorithm>
#include <iostream>

#include "abpoa.h"

// FIXME: we already have something similar in PingPong.hpp - merge?
// AaCcGgTtNn ==> 0,1,2,3,4
extern const unsigned char _nt4_table[256];

struct Fragment {
    std::string name;
    uint ref_s;
    uint ref_e;
    uint read_s;
    uint read_e;
    uint o_ref_s ;
    uint o_ref_e ;
    bool left_extended = false ;
    bool right_extended = false ;
    uint t; // 0: only D; 1: only I; 2: D and I
    std::string seq ;
    std::string qual ;
    std::vector<std::pair<int, int>> aligned_pairs ;

    Fragment() {}

    Fragment(const std::string& name_, uint _ref_s, uint _ref_e, uint _read_s, uint _read_e, uint t_, const std::string& seq_ = "", const std::string& qual_ = "") {
        name = name_ ;
        ref_s = _ref_s ;
        ref_e = _ref_e ;
        //
        o_ref_s = _ref_s ;
        o_ref_e = _ref_e ;
        //
        read_s = _read_s ;
        read_e = _read_e ;
        t = t_ ;
        seq = seq_ ;
        qual = qual_ ;
    }

    uint size() const {
        return seq.size();
    }

    bool operator<(const Fragment &f) const {
        return ref_s < f.ref_s;
    }

    bool operator==(const Fragment &f) const {
        return name == f.name and ref_s == f.ref_s and ref_e == f.ref_e;
    }

    void extend(bool left, bool right) {
        left_extended = left ;
        right_extended = right ;
    }
};

class Cluster {
public:
    std::vector<Fragment> fragments;
    std::string chrom;
    uint s;
    uint e;
    uint full_cov ;

    Cluster(const std::string &);
    void add_fragment(Fragment);
    void set_full_coverage(const uint);
    std::string poa() const;
    uint get_type() const;
    void clear();
    uint size() const;
    bool empty() const;
    std::string get_id() const ;
    const Fragment &front() const;
    const Fragment &back() const;
    std::vector<Fragment>::const_iterator begin() const;
    std::vector<Fragment>::const_iterator end() const;

    bool operator<(const Cluster &c) const {
        return s < c.s;
    }

    bool operator==(const Cluster &c) const {
        return chrom == c.chrom and s == c.s and e == c.e ;
    }
};

#endif
