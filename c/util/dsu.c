#include <stdlib.h>
#include <stddef.h>
#include <stdint.h>
#include "dsu.h"

void dsu_init(dsu* u, uint64_t n) {
  u->e = (int64_t*) malloc(sizeof(int64_t) * n);
  u->n = n;
  for(size_t i = 0; i < n; i++) u->e[i] = i;
}

uint64_t dsu_find_set(dsu* u, uint64_t x) {
  if(u->e[x] == x) return x;
  u->e[x] = dsu_find_set(u, u->e[x]);
  return u->e[x];
}

int8_t dsu_union_set(dsu* u , uint64_t x, uint64_t y) {
  x = dsu_find_set(u, x);
  y = dsu_find_set(u, y);
  if(x == y) return 0;
  if(x > y) {
    int64_t aux = x;
    x = y;
    y = aux;
  }
  u->e[y] = x;
  return 1;
}

void dsu_free(dsu* u) {
  free(u->e);
}

