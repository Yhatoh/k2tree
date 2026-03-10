#ifndef __DSU_H__
#define __DSU_H__

#include <stdlib.h>
#include <stddef.h>
#include <stdint.h>

typedef struct {
  uint64_t n;
  uint64_t* e;
} dsu;

void dsu_init(dsu* u, uint64_t n);
uint64_t dsu_find_set(dsu* u, uint64_t x);
int8_t dsu_union_set(dsu* u , uint64_t x, uint64_t y);
void dsu_free(dsu* u);

#endif
