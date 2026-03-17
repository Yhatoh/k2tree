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

void bv_init(bv_t* z) {
  z->n = 0;
  z->maxn = 5; // 5 * 64 elements
  z->a = (uint64_t*) malloc(sizeof(uint64_t) * z->maxn);
}

void bv_pb(bv_t* z, uint64_t bit) {
  if(z->n >= z->maxn * 64) {
    z->maxn *= 2;
    z->a = (uint64_t*) realloc(z->a, sizeof(uint64_t) * z->maxn);
    if(z->a == NULL) quit("realloc failed",__LINE__,__FILE__);
  }


  z->a[z->n / 64] = (bit << (z->n % 64)) | (z->a[z->n / 64] & ((1LL << (z->n % 64)) - 1));
  z->n++;
}

void bv_append_int(bv_t* z, uint64_t num) {
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
                    //(z->a[z->n / 64] & ((1LL << (64 - z->n % 64)) - 1));
                    (z->a[z->n / 64] & ((1LL << (z->n % 64)) - 1));
  z->a[z->n / 64 + 1] = num >> (64 - z->n % 64);
  z->n += 64;
}

void bv_append_uw(bv_t* z, uint64_t num, uint64_t w) {
  bv_append_int(z, num);
  z->n -= 64;
  z->n += w;
}

void bv_append_u16(bv_t* z, uint64_t num) {
  if(z->n + 16 >= z->maxn * 64) {
    z->maxn *= 2;
    z->a = (uint64_t*) realloc(z->a, sizeof(uint64_t) * z->maxn);
    if(z->a == NULL) quit("realloc failed",__LINE__,__FILE__);
  }
  bv_append_int(z, num);
  z->n -= 64;
  z->n += 16;
}

void bv_append(bv_t* a, bv_t* b) {
  size_t i = 0;
  for(; i + 64 < b->n; i += 64) {
    bv_append_int(a, bv_get_int(b, i, 64));
  }

  for(; i < b->n; i++) {
    bv_pb(a, bv_i(b, i));
  }
}


void bv_free(bv_t* z) {
  if(z->a != NULL)
    free(z->a);
  z->a = NULL;
  z->n = 0;
  z->maxn = 0;
}

void bv_shrink(bv_t* z) {
  z->maxn = (z->n + 64 - 1) / 64;
  z->a = (uint64_t*) realloc(z->a, sizeof(uint64_t) * z->maxn);
  if(z->a == NULL) quit("realloc failed",__LINE__,__FILE__);
}

uint64_t bv_i(const bv_t* z, size_t i) {
  assert(i < z->n);
  return (z->a[i / 64] & (1LL << (i % 64))) >> (i % 64);
}

uint64_t bv_get_int(const bv_t* z, size_t i, uint8_t len) {
  assert(len >= 1 && len <= 64);
  if((i % 64) + len > 64) {
    // divided in two
    return (z->a[i / 64] >> (i % 64)) |
      ((z->a[(i / 64) + 1] & ((1LL << (len - (64 - (i % 64)))) - 1)) << (64 - (i % 64)));
  }
  if(len == 64) return z->a[i / 64];
  return (z->a[i / 64] >> (i % 64)) & ((1LL << len) - 1);
}

void bv_reserve(bv_t* z, size_t m) {
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

void bv_grow(bv_t* z, size_t i) {
  z->n += i;
  while(z->n > z->maxn * 64) {
    z->maxn *= 2;
    z->a = (uint64_t*) realloc(z->a, sizeof(uint64_t) * z->maxn);
    if(z->a == NULL) quit("realloc failed",__LINE__,__FILE__);
  }
}

uint64_t bv_size_in_bits(const bv_t* z) {
  return sizeof(uint64_t) * z->maxn + sizeof(bv_t);
}

void bv_save_to_file(const bv_t* z, const char* fname) {
  FILE* out = fopen(fname, "w");
  fwrite(&(z->n), sizeof(size_t), 1, out);
  size_t w = fwrite(z->a, sizeof(uint64_t), (z->n + 64 - 1) / 64, out);
  if(w != (z->n + 64 - 1) / 64) quit("bv_save_to_file: error writing bv_t to file",__LINE__,__FILE__);
  fclose(out);
}

void bv_load_from_file(bv_t* z, const char* fname) {
  FILE* in = fopen(fname, "r");
  size_t w = fread(&(z->n), sizeof(size_t), 1, in);
  if(w != 1) quit("bv_load_from_file: error reading n from file", __LINE__, __FILE__);
  z->maxn = (z->n + 64 - 1) / 64;
  z->a = (uint64_t*) malloc(sizeof(uint64_t) * z->maxn);
  w = fread(z->a, sizeof(size_t), z->maxn, in);
  if(w != z->maxn) quit("bv_load_from_file: error reading n from file", __LINE__, __FILE__);
  fclose(in);
}

// write error message and exit
static void quit(const char *msg, int line, char *file) {
  if(errno==0)  fprintf(stderr,"== %d == %s\n",getpid(), msg);
  else fprintf(stderr,"== %d == %s: %s\n",getpid(), msg,
               strerror(errno));
  fprintf(stderr,"== %d == Line: %d, File: %s\n",getpid(),line,file);
  exit(1);
}
