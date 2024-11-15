#ifndef UTIL_HPP
#define UTIL_HPP
#include<cstdint>

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
