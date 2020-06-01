#ifndef _BLOOM_FILTER_HPP
#define _BLOOM_FILTER_HPP

#include <bitset>
#include <vector>
#include <time.h>

#include "xxhash.hpp"

using namespace std;

template <uint32_t size, uint8_t nhashes>
class BF {
public:
  BF() : _seeds(nhashes, 0) {
    srand48(time(NULL));
    _size = size;
    for(uint i=0; i<_seeds.size(); ++i)
      _seeds[i] = lrand48();
  }
  ~BF() {}

  BF(const BF &bf) : _seeds(bf._seeds.size(), 0) {
    _size = bf._size;
    for(uint i=0; i<_seeds.size(); ++i)
      _seeds[i] = bf._seeds[i];
  }

  void add(const uint64_t &S) {
    for(const uint seed : _seeds) {
      uint64_t h = xxh::xxhash<64>(&S, sizeof(uint64_t), seed);
      _bv.set(h%_size);
    }
  }
  
  bool test(const uint64_t S) const {
    bool res = 1;
    for(const uint seed : _seeds) {
      uint64_t h = xxh::xxhash<64>(&S, sizeof(uint64_t), seed);
      res &= _bv[h%_size];
    }
    return res;
  }


private:
  const BF &operator=(const BF &other) { return *this; }
  const BF &operator=(const BF &&other) { return *this; }

  bitset<size> _bv;
  uint32_t _size;
  vector<int32_t> _seeds;
};

#endif
