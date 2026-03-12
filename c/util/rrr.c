#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include "rrr.h"
#include "../k2bp/util.h"
#include "rrr_util.h"

void rrr_compress(rrr_t* rrr, uint8_t b, uint8_t* bits, size_t n) {
  rrr->n = n;
  rrr->b = b;
  rrr->logb = ceil_log2(b + 1);
  size_t pairs = (n + b - 1) / b;
  size_t pairs_per_64 = 64 / rrr->logb;
  iv_init(&(rrr->c), pairs, rrr->logb);
  bv_init(&(rrr->o));
  size_t curr_bit_o = 0;
  size_t curr_pair = 0;
  for(size_t i = 0; i < n; i += 63) {
    size_t class = 0;
    for(size_t j = i; j < i + 63 && j < n; j++) {
      if(bits[j / 8] & (1 << (j % 8))) 
        class++;
    }

    size_t class_ = class;
    size_t offset = 0;
    size_t j = 0;
    while(class_ > 0 && class_ <= b - j) {
      if(bits[(j + i) / 8] & (1 << ((j + i) % 8))) {
        offset += comb[b - j - 1][class_];
        class_--;
      }
      j++;
    }
    iv_set(&(rrr->c), curr_pair, class);
    bv_append_uw(&(rrr->o), offset, l[class]);
    curr_pair++;
  }
  bv_shrink(&(rrr->o));
}

void rrr_free(rrr_t* rrr) {
  rrr->n = rrr->b = 0;
  iv_free(&(rrr->c));
  bv_free(&(rrr->o));
}

void rrr_decompress(rrr_t* rrr, uint8_t* bits, size_t* n) {
}
