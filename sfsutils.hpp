#ifndef SFSUTILS_HPP
#define SFSUTILS_HPP

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

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
  uint s;
  uint l;
  uint c;

  SFS() {
    s = 0;
    l = 0;
    c = 0;
  }

  SFS(uint s_, uint l_, uint c_) {
    s = s_;
    l = l_;
    c = c_;
  }

  void reverse(uint p) { s = p - s - l; }
};

bool operator<(const SFS &, const SFS &);

map<string, vector<SFS>> parse_sfsfile(const string &, int);

#endif
