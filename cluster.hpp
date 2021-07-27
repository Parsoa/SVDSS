#ifndef CLUSTER_HPP
#define CLUSTER_HPP

#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <algorithm>

#include "abpoa.h"

using namespace std;

// FIXME: we already have something similar in PingPong.hpp - merge?
// AaCcGgTtNn ==> 0,1,2,3,4
extern const unsigned char _nt4_table[256];

// CHECKME: maybe we already have something like this (a fragment is a superstring)
struct Fragment
{
    string name;
    uint s;
    uint e;
    uint t; // 0: only D; 1: only I; 2: D and I
    string seq;

    Fragment() {}

    Fragment(const string &name_, uint s_, uint e_, uint t_, const string &seq_ = "")
    {
        name = name_;
        s = s_;
        e = e_;
        t = t_;
        seq = seq_;
    }

    uint size() const
    {
        return seq.size();
    }

    bool operator<(const Fragment &f) const
    {
        return s < f.s;
    }

    bool operator==(const Fragment &f) const
    {
        return name == f.name and s == f.s and e == f.e;
    }
};

class Cluster
{
public:
    list<Fragment> fragments;
    string chrom;
    uint s;
    uint e;
    uint full_cov;

    Cluster(const string &);
    void add_fragment(const Fragment &);
    void set_full_coverage(const uint);
    string poa() const;
    uint get_type() const;
    void clear();
    uint size() const;
    bool empty() const;
    const Fragment &front() const;
    const Fragment &back() const;
    list<Fragment>::const_iterator begin() const;
    list<Fragment>::const_iterator end() const;
};

#endif