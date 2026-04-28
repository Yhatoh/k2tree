#include <assert.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include "iv.h"

void iv_init(iv_t* z, size_t n, uint8_t w) {
  assert(w > 0 && w <= 64);

  z->n = n;
  z->w = w;
  z->data = (uint64_t*) malloc(sizeof(uint64_t) * (((n * w) + 64 - 1) / 64 + 2));
}

uint64_t iv_get(const iv_t* z, size_t i) {
  assert(i >= 0 && i < z->n);

  size_t bit = i * z->w;
  size_t block = bit / 64;
  size_t offset = bit % 64;

  size_t ret = z->data[block] >> offset;
  if(offset + z->w > 64) {
    ret |= z->data[block + 1] << (64 - offset);
  }

  return ret & ((1ULL << z->w) - 1);
}

void iv_set(iv_t* z, size_t i, uint64_t num) {
  assert(i >= 0 && i < z->n);

  size_t bit = i * z->w;
  size_t block = bit / 64;
  size_t offset = bit % 64;

  z->data[block] = (z->data[block] & ((1ULL << offset) - 1)) | (num << offset);
  if(offset + z->w > 64) {
    z->data[block + 1] = (num >> (64 - offset)) & ((1ULL << (z->w - 64 - offset)) - 1);
  }
}

void iv_free(iv_t* z) {
  z->n = 0;
  z->w = 0;
  free(z->data);
}

void iv_copy(const iv_t* src, iv_t* dest) {
  dest->n = src->n;
  dest->w = src->w;
  dest->data = (uint64_t*) malloc(sizeof(uint64_t) * (((src->n*src->w)+64-1)/64+2));
  for(size_t i = 0; i < (((src->n*src->w)+64-1)/64+2); i++)
    dest->data[i] = src->data[i];
}
