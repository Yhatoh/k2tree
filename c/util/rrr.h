#ifndef __RRR_H__
#define __RRR_H__

#include <stddef.h>
#include <stdint.h>
#include "iv.h"
#include "bv_t.h"

typedef struct rrr_t {
  size_t n;
  uint8_t b;
  uint64_t logb;
  iv_t c;
  bv_t o;
} rrr_t;

void rrr_compress(rrr_t* rrr, uint8_t b, uint8_t* bits, size_t n);
void rrr_free(rrr_t* rrr);
void rrr_decompress(rrr_t* rrr, uint8_t* bits, size_t* n);

#endif
