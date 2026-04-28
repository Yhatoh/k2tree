#ifndef __IV_H__
#define __IV_H__

#include <stddef.h>
#include <stdint.h>

typedef struct iv_t {
  size_t n;
  uint8_t w;
  uint64_t *data;
} iv_t;

void iv_init(iv_t* z, size_t n, uint8_t w);
uint64_t iv_get(const iv_t* z, size_t i);
void iv_set(iv_t* z, size_t i, uint64_t num);
void iv_free(iv_t* z);
void iv_copy(const iv_t* src, iv_t* dest);

#endif
