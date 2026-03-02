#include <assert.h>
#include <stdlib.h>
#include <errno.h>
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <string.h>
#include <unistd.h>

#include "bv_t.h"

static void quit(const char *msg, int line, char *file);

void print_bits(uint64_t num) {
  for(size_t i = 0; i < 64; i++) {
    if(num & 1) printf("1");
    else printf("0");
    num >>= 1;
  }
  printf("\n");
}

void bv_t_init(bv_t* z) {
  z->n = 0;
  z->maxn = 5; // 5 * 64 elements
  z->a = (uint64_t*) malloc(sizeof(uint64_t) * z->maxn);
}

void bv_t_pb(bv_t* z, uint64_t bit) {
  if(z->n >= z->maxn * 64) {
    z->maxn *= 2;
    z->a = (uint64_t*) realloc(z->a, sizeof(uint64_t) * z->maxn);
    if(z->a == NULL) quit("realloc failed",__LINE__,__FILE__);
  }


  z->a[z->n / 64] = (bit << (z->n % 64)) | (z->a[z->n / 64] & ((1LL << (z->n % 64)) - 1));
  z->n++;
}

void bv_t_append_int(bv_t* z, uint64_t num) {
  if(z->n + 64 >= z->maxn * 64) {
    z->maxn *= 2;
    z->a = (uint64_t*) realloc(z->a, sizeof(uint64_t) * z->maxn);
    if(z->a == NULL) quit("realloc failed",__LINE__,__FILE__);
  }
  if(z->n % 64 == 0) {
    z->a[z->n / 64] = num;
    z->n += 64;
    return;
  }

  z->a[z->n / 64] = (num << z->n % 64) |
                    (z->a[z->n / 64] & ((1LL << (64 - z->n % 64)) - 1));
  z->a[z->n / 64 + 1] = num >> (64 - z->n % 64);
  z->n += 64;
}

void bv_t_append(bv_t* a, bv_t* b) {
  size_t i = 0;
  for(; i < b->n / 64; i++) {
    bv_t_append_int(a, b->a[i]);
  }

  if(b->n % 64 > 0) {
    bv_t_append_int(a, b->a[b->n / 64]);
    a->n -= 64 - b->n % 64;
  }
}


void bv_t_free(bv_t* z) {
  free(z->a); z->a = NULL;
  z->n = 0;
  z->maxn = 0;
}

void bv_t_shrink(bv_t* z) {
  z->maxn = (z->n + 64 - 1) / 64;
  z->a = (uint64_t*) realloc(z->a, sizeof(uint64_t) * z->maxn);
  if(z->a == NULL) quit("realloc failed",__LINE__,__FILE__);
}

uint64_t bv_t_i(bv_t* z, size_t i) {
  assert(i < z->n);
  return (z->a[i / 64] & (1LL << (i % 64))) >> (i % 64);
}

uint64_t bv_t_get_int(bv_t* z, size_t i, uint8_t len) {
  assert(len >= 1 && len <= 64);
  if((i % 64) + len > 64) {
    // divided in two
    return (z->a[i / 64] >> (i % 64)) |
      ((z->a[(i / 64) + 1] & ((1LL << (len - (64 - (i % 64)))) - 1)) << (64 - (i % 64)));
  }
  if(len == 64) return z->a[i / 64];
  return (z->a[i / 64] >> (i % 64)) & ((1 << len) - 1);
}

void bv_t_reserve(bv_t* z, size_t m) {
  if(z->a != NULL) {
    z->a = (uint64_t*) realloc(z->a, sizeof(uint64_t) * m);
    if(z->a == NULL) quit("realloc failed",__LINE__,__FILE__);
    z->maxn = m;
    if(z->n > z->maxn) z->n = z->maxn * 64; // erase values that doesn't exist anymore
    return;
  }

  z->a = (uint64_t*) malloc(sizeof(uint64_t) * m);
  if(z->a == NULL) quit("malloc failed",__LINE__,__FILE__);
  z->maxn = m;
}

void bv_t_grow(bv_t* z, size_t i) {
  z->n += i;
  while(z->n > z->maxn * 64) {
    z->maxn *= 2;
    z->a = (uint64_t*) realloc(z->a, sizeof(uint64_t) * z->maxn);
    if(z->a == NULL) quit("realloc failed",__LINE__,__FILE__);
  }
}

uint64_t size_in_bits(bv_t* z) {
  return sizeof(uint64_t) * z->maxn + sizeof(bv_t);
}

// write error message and exit
static void quit(const char *msg, int line, char *file) {
  if(errno==0)  fprintf(stderr,"== %d == %s\n",getpid(), msg);
  else fprintf(stderr,"== %d == %s: %s\n",getpid(), msg,
               strerror(errno));
  fprintf(stderr,"== %d == Line: %d, File: %s\n",getpid(),line,file);
  exit(1);
}
