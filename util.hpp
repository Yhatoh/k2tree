#ifndef UTIL_HPP
#define UTIL_HPP

#include <cstdint>
#include <queue> 
#include <vector>
#include <utility>       
#include <map>
#include <cmath>
#include <inttypes.h>
#include <string>
#include <tuple>

#include <sdsl/vectors.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/coder_elias_delta.hpp>
#include <sdsl/suffix_arrays.hpp>

using namespace std;
// from succint repository
static const uint8_t debruijn64_mapping[64] = {
  63,  0, 58,  1, 59, 47, 53,  2,
  60, 39, 48, 27, 54, 33, 42,  3,
  61, 51, 37, 40, 49, 18, 28, 20,
  55, 30, 34, 11, 43, 14, 22,  4,
  62, 57, 46, 52, 38, 26, 32, 41,
  50, 36, 17, 19, 29, 10, 13, 21,
  56, 45, 25, 31, 35, 16,  9, 12,
  44, 24, 15,  8, 23,  7,  6,  5
};

static const uint64_t debruijn64 = 0x07EDD5E59A4E28C2ULL;

inline uint8_t bit_position(uint64_t x){
    return debruijn64_mapping[(x * debruijn64) >> 58];
}

inline uint8_t msb(uint64_t x, unsigned long& ret){
  if (!x)
    return false;

  x |= x >> 1;
  x |= x >> 2;
  x |= x >> 4;
  x |= x >> 8;
  x |= x >> 16;
  x |= x >> 32;

  x ^= x >> 1;
  ret = bit_position(x);

  return true;
}

static const uint8_t table_mul[16][16] = {{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
  {0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3},
  {0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3},
  {0,1,2,3,1,1,3,3,2,3,2,3,3,3,3,3},
  {0,4,8,12,0,4,8,12,0,4,8,12,0,4,8,12},
  {0,5,10,15,0,5,10,15,0,5,10,15,0,5,10,15},
  {0,4,8,12,1,5,9,13,2,6,10,14,3,7,11,15},
  {0,5,10,15,1,5,11,15,2,7,10,15,3,7,11,15},
  {0,0,0,0,4,4,4,4,8,8,8,8,12,12,12,12},
  {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15},
  {0,0,0,0,5,5,5,5,10,10,10,10,15,15,15,15},
  {0,1,2,3,5,5,7,7,10,11,10,11,15,15,15,15},
  {0,4,8,12,4,4,12,12,8,12,8,12,12,12,12,12},
  {0,5,10,15,4,5,14,15,8,13,10,15,12,13,14,15},
  {0,4,8,12,5,5,13,13,10,14,10,14,15,15,15,15},
  {0,5,10,15,5,5,15,15,10,15,10,15,15,15,15,15}
};

inline uint8_t minimat_mul(uint8_t a, uint8_t b) { return table_mul[a][b]; }

static const uint8_t rev_table[16] = {0,8,4,12,2,10,6,14,1,9,5,13,3,11,7,15};

inline uint8_t msb(uint64_t x){
  unsigned long ret = -1U;
  msb(x, ret);
  return (uint8_t)ret;
}

inline uint64_t ceil_log2(uint64_t x) {
  return (x > 1) ? msb(x - 1) + 1 : 0;
}

inline uint64_t floor_log2(uint64_t x) {
  return (x > 1) ? msb(x) : 0;
}

#endif
