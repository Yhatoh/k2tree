#ifndef __BV_H__
#define __BV_H__

#include <stdint.h>
#include <stddef.h>

typedef struct bv_t {
  size_t n;
  size_t maxn;
  uint64_t* a;
} bv_t;

void bv_t_init(bv_t* z);
void bv_t_pb(bv_t* z, uint64_t bit); // tested
void bv_t_append(bv_t* a, bv_t* b); // tested
void bv_t_append_int(bv_t* z, uint64_t num); // tested
void bv_t_free(bv_t* z); // should be correct
void bv_t_reserve(bv_t* z, size_t m); // should be correct
void bv_t_grow(bv_t* z, size_t i); // should be correct
void bv_shrink(bv_t* z); // should be correct
uint64_t bv_t_i(bv_t* z, size_t i); // tested
uint64_t bv_t_get_int(bv_t* z, size_t i, uint8_t len); // tested
uint64_t size_in_bits(bv_t* z);


#endif
