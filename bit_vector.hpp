#ifndef __BIT_VECTOR__
#define __BIT_VECTOR__

// c++ includes
#include <cassert>
#include <cstdint>
#include <vector>


struct bvector {
  std::vector< uint8_t > bv;
  uint64_t n;

  bvector() {
    n = 0;
  }

  uint64_t size() {
    return n;
  }

  uint8_t operator[](const uint64_t i) {
    assert(i <= n);
    const uint64_t integer = i / 8;
    const uint64_t i_int = i % 8;

    return ((bv[integer] & (1 << i_int)) != 0);
  }

  void set(const uint64_t num, const uint64_t i, const uint64_t w) {
    assert(i + w <= n);
    uint64_t integer = i / 8;
    uint64_t i_int = i % 8;

    for(uint64_t bit = 0; bit < w; bit++) {
      if(num & (1 << bit)) {
        bv[integer] |= 1 << i_int;
      } else {
        bv[integer] &= ~(1 << i_int);
      }
      i_int++;
      if(i_int == 8) {
        integer++;
        i_int = 0;
      }
    }
  }
    
  void reserve(uint64_t i) {
    bv.reserve((i + 8 - 1) / 8);
  }

  void push_back(uint8_t bit) {
    if(n / 8 >= bv.size()) {//full
      bv.push_back(0);
    }
    const uint64_t integer = n / 8;
    const uint64_t i_int = n % 8;
    n++;
    bv[integer] |= (1 << i_int);
  }

  void destroy() {
    std::vector< uint8_t >().swap(bv);
  }

  void concat(bvector& a, uint64_t i = 0, uint64_t end = 0) {
    if(end == 0) end = a.size();
    while(i < end && 
          !(i % 8 == 0 && size() % 8 == 0)) {
      push_back(a.bv[i]);
      i++;
    }

    if(i == a.size()) return;

    auto itr_begin = a.bv.begin() + (i / 8);
    bv.insert(bv.end(), itr_begin, a.bv.end());
    n += a.size() - i / 8;
  }

  void clear() {
    bv.clear();
    n = 0;
  }
};

#endif // !__BIT_VECTOR__

