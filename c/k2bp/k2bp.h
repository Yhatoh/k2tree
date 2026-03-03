#ifndef __K2BP_H__
#define __K2BP_H__

#include <stdint.h>
#include "bv_t.h"

#define NUM_SUPPORT 36
#define GET_NODES(x) (x & ((1LL << NUM_SUPPORT) - 1))
#define GET_SKIPS(x) (x >> NUM_SUPPORT)
#define ENCODE(x, y) (x << NUM_SUPPORT) | y

#define _K_ 2
#define BLOCK_SIZE 256
#define SAMPLE_SIZE(n) ((n + BLOCK_SIZE - 1) / BLOCK_SIZE + 1)

// usefull information during traversals
typedef struct k2bp_traversal_t {
  size_t msize;
  size_t x, y;
  size_t i_t;
  size_t i_l;
} k2bp_traversal_t;

typedef struct k2bp_t {
  size_t msize; // pow 2 matrix size
  size_t rmsize; // real matrix size
  size_t m; // amount of ones

  // tree
  bv_t t;

  // array l
  size_t maxn_l;
  size_t n_l;
  uint8_t* l;

  // tree support helper
  //  child support
  //  amount of leaves for a subtree
  uint16_t* exc_min_samples;
  uint16_t* exc_samples;
  uint8_t* leaves_samples; // this is specific for block 256

  // subtree information
  size_t n_info;
  size_t threshold;
  uint64_t* subtreeinfo;
  uint32_t* leavesinfo;
} k2bp_t;

#define K2BP_INITIALIZER {0, 0, 0, {0, 0, NULL}, 0, 0, NULL, NULL, NULL, NULL, 0, 0, NULL, NULL}

// k2 tree operations
size_t k2bp_build_from_textfile(k2bp_t* a, const char* f, size_t fsize);
void k2bp_write_leaf(k2bp_t* a, uint8_t leaf);
uint8_t k2bp_read_leaf(const k2bp_t*a, const size_t pos);
size_t k2bp_compute_height(const size_t rmsize);
void k2bp_free(k2bp_t* a);
void k2bp_dfs(k2bp_traversal_t* pos_a, const k2bp_t* a, size_t* nodes, size_t* leaves, size_t* nz, size_t* levels, size_t curr_level);

size_t k2bp_show_stats(const k2bp_t* a, const char* fname, FILE* f);
size_t k2bp_stats(const k2bp_t* a, size_t* nodes, size_t* leaves, size_t* nz);

void k2bp_save_to_file(const k2bp_t* a, const char* fname);
void k2bp_load_from_file(k2bp_t* a, const char* fname);

// matrix operations
uint32_t* k2bp_nonzeros(const k2bp_t* a, size_t* n);

#endif
