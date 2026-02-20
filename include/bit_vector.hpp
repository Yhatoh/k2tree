#ifndef __BIT_VECTOR__
#define __BIT_VECTOR__

// c++ includes
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "debug.hpp"
// super simple append only bit_vector

const size_t b64 = 6;
const size_t m64 = 63;

struct bvector {
  std::vector< uint64_t > bv;
  uint64_t n;

  bvector() {
    n = 0;
  }

  uint64_t size() {
    return n;
  }

  uint64_t operator[](const uint64_t i) {
    assert(i <= n);
    const uint64_t integer = i >> b64;
    const uint64_t i_int = i & m64;

    return ((bv[integer] & (1LL << i_int)) > 0);
  }

  void set(const uint64_t num, const uint64_t i, const uint64_t w) {
    assert(i + w <= n);
    uint64_t integer = i >> b64;
    uint64_t i_int = i & m64;

    for(uint64_t bit = 0; bit < w; bit++) {
      if(num & (1LL << bit)) {
        bv[integer] |= 1LL << i_int;
      } else {
        bv[integer] &= ~(1LL << i_int);
      }
      i_int++;
      if(i_int == 8) {
        integer++;
        i_int = 0;
      }
    }
  }
    
  void reserve(uint64_t i) {
    bv.reserve((i + 64 - 1) >> b64);
  }

  void push_back(uint64_t bit) {
    if((n & m64) == 0) {//full
      bv.push_back(0);
    }
    const uint64_t integer = n >> b64;
    const uint64_t i_int = n & m64;
    n++;
    bv[integer] |= (bit << i_int);
  }

  void destroy() {
    std::vector< uint64_t >().swap(bv);
  }
  
  uint64_t read_int(uint64_t i, uint64_t len = 64) {
    assert(len > 0);

    if(((i & m64) + len) > 64) {
      // has to be divided in two
      return (bv[i >> b64] >> (i & m64)) |
             ((bv[(i >> b64) + 1] & ((1LL << (len - (64 - (i & m64)))) - 1)) << (64 - (i & m64)));
    }
    if(len == 64) return bv[i >> b64];
    return (bv[i >> b64] >> (i & m64)) & ((1LL << len) - 1);
  }

  // this functions is super slow
  void concat(bvector& a, uint64_t i = 0, uint64_t end = 0) {
    assert(i <= a.size());
    assert(end <= a.size());
    if(end == 0) end = a.size();

    for(; (n & m64) != 0 && (i < end); i++)
      push_back(a[i]);

    for(; i + 64 < end; i += 64) {
      bv.push_back(a.read_int(i, 64));
      n += 64;
    }

    if(end - i > 0) {
      bv.push_back(a.read_int(i, end - i));
      n += end - i;
    }
  }

  void clear() {
    bv.clear();
    n = 0;
  }
};

#endif // !__BIT_VECTOR__

