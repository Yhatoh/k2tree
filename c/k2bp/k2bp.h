#ifndef __K2BP_H__
#define __K2BP_H__

#include <stdint.h>
#include "../util/bv_t.h"
#include "../util/iv.h"
#include "../util/rrr.h"

#define NUM_SUPPORT 32
#define GET_NODES(x) (x & ((1ULL << NUM_SUPPORT) - 1))
#define GET_SKIPS(x) (x >> NUM_SUPPORT)
#define ENCODE(x, y) (x << NUM_SUPPORT) | y
#define HAS_POINTERS(a) (a->pointers.data != NULL)

#define _K_ 2
#define BLOCK_SIZE 256
#define SAMPLE_SIZE(n) ((n + BLOCK_SIZE - 1) / BLOCK_SIZE + 1)

#define SUPER_BLOCK_SIZE 8192 // extra fast traversal, specific for compressed BP
#define SAMPLE_SIZE_SUPER(n) ((n + SUPER_BLOCK_SIZE - 1) / SUPER_BLOCK_SIZE + 1)

#define BLOCK_SIZE_RANK 4096 // apparently this make slow 
                             // faster matrix-matrix multiplication, this should be reduce
#define SAMPLE_SIZE_RANK(n) ((n + BLOCK_SIZE_RANK - 1) / BLOCK_SIZE_RANK + 1)

#define LEAF_0 1
#define LEAF_1 3
#define LEAF_P 11

// usefull information during traversals
typedef struct k2bp_traversal_t {
  size_t msize;
  size_t x, y;
  size_t i_t;
  size_t i_l;
  size_t node;
  uint64_t size;
  uint32_t leaves;
  int16_t excess;
  size_t i_p;
  uint8_t flag_p;
//  uint32_t* rank_pointers;
//  uint32_t* leaves_pointers;
//  uint32_t* node_pointers;
//  uint64_t* size_sub_pointers;
  iv_t rank_pointers;
  iv_t leaves_pointers;
//  iv_t node_pointers;
  uint8_t flag_cl;
} k2bp_traversal_t;

//#define K2BP_TRAVERSAL_INITIALIZER {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, NULL, NULL, NULL, NULL, 0}
//#define K2BP_TRAVERSAL_INITIALIZER {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, {0, 0, NULL}, {0, 0, NULL}, {0, 0, NULL}, 0}
#define K2BP_TRAVERSAL_INITIALIZER {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, {0, 0, NULL}, {0, 0, NULL}, 0}

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

  rrr_t cl;

  // tree support helper
  //  child support
  //  amount of leaves for a subtree
  uint16_t* exc_min_samples; // I think this two can be reduce to uint8_t test it later
  uint16_t* exc_samples;
  uint8_t* leaves_samples; // this is specific for block <= 256
                           // if you want bigger blocks you should change this type
  uint8_t* pointers_samples; // only for compressed version

  // subtree information
  size_t n_info;
  size_t threshold;
  uint64_t* subtreeinfo;
  uint32_t* leavesinfo;

  // pointers info
  size_t n_p;
  iv_t pointers;

  // extra sampling for compressed version
  uint16_t* super_exc_min_samples; // I think this two can be reduce to uint8_t test it later
  uint16_t* super_exc_samples;
  uint16_t* super_leaves_samples;
  uint16_t* super_pointers_samples;

} k2bp_t;

#define K2BP_INITIALIZER {0, 0, 0, {0, 0, NULL}, 0, 0, NULL, \
                          {0, 0, 0, {0, 0, NULL}, {0, 0, NULL}}, \
                          NULL, NULL, NULL, NULL, 0, 0, NULL, NULL, \
                          0, {0, 0, NULL}, NULL, NULL, NULL, NULL}

// k2 tree operations
size_t k2bp_build_from_textfile(k2bp_t* a, const char* f, size_t fsize);
void k2bp_write_leaf(k2bp_t* a, uint8_t leaf);
uint8_t k2bp_read_leaf(const k2bp_t*a, const size_t pos);
size_t k2bp_compute_height(const size_t rmsize);
void k2bp_free(k2bp_t* a);
void k2bp_dfs(k2bp_traversal_t* pos_a, const k2bp_t* a,
              size_t* nodes, size_t* leaves, size_t* nz, size_t* levels, size_t curr_level);
void k2bp_build_exc_sampling(k2bp_t* a);
void k2bp_addsubtree_info(k2bp_t* a, size_t threshold);
size_t k2bp_checksubtree_info(k2bp_t* a);
uint8_t k2bp_equal(k2bp_t* a, k2bp_t* b);
void k2bp_compress_subtrees(k2bp_t* a, k2bp_t* c, size_t limit);
void k2bp_decompress_subtrees(k2bp_t* c, k2bp_t* a);
void k2bp_compress_leaves(k2bp_t* a);
void k2bp_decompress_leaves(k2bp_t* a);
void k2bp_copy(k2bp_traversal_t* pos_src, const k2bp_t* src, k2bp_t* dest);

// k2 tree information
size_t k2bp_show_stats(k2bp_t* a, const char* fname, FILE* f);
size_t k2bp_stats(k2bp_t* a, size_t* nodes, size_t* leaves, size_t* nz);

// io k2 tree
void k2bp_save_to_file(const k2bp_t* a, const char* fname);
void k2bp_load_from_file(k2bp_t* a, const char* fname);

// matrix operations
uint32_t* k2bp_nonzeros(k2bp_t* a, size_t* n);
void k2bp_sum(k2bp_t* a, k2bp_t* b, k2bp_t* c);
void k2bp_scansum(k2bp_t* a, k2bp_t* b, k2bp_t* c);
void k2bp_mul(k2bp_t* a, k2bp_t* b, k2bp_t* c);
void k2bp_scanmul(k2bp_t* a, k2bp_t* b, k2bp_t* c);

#endif
