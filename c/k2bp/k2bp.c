#include <assert.h>
#include <limits.h>
#include <errno.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "k2bp.h"
#include "util.h"
#include "../util/bv_t.h"
#include "../util/iv.h"
#include "../util/vu64.h"
#include "../util/dsu.h"
#include "../../libsais/include/libsais64.h"

static size_t binsearch(uint64_t *ia, size_t n, uint64_t x);
static uint64_t *create_ia(FILE *f, size_t *n, size_t *msize, size_t xsize);
static size_t create_k2bp(uint64_t ia[], size_t n, size_t msize, k2bp_t *a);
static void quit(const char *msg, int line, char *file);
static uint8_t encode_leaf(uint64_t ia[], size_t n, size_t smin);
static void reck2bp_nonzeros(k2bp_traversal_t* pos_a, const k2bp_t* a, uint32_t* arr, size_t* n);
static void reck2bp_addsubtree_info(k2bp_traversal_t* pos_a, k2bp_t* a, uint32_t* leaves, vu64_t* subinfo, vu64_t *leavesinfo);
static void reck2bp_checksubtree_info(k2bp_traversal_t* pos_a, const k2bp_t* a);

static void reck2bp_sum(k2bp_traversal_t* pos_a, const k2bp_t* a, k2bp_traversal_t* pos_b, const k2bp_t* b, k2bp_t* c);
static void reck2bp_scansum(k2bp_traversal_t* pos_a, const k2bp_t* a, k2bp_traversal_t* pos_b, const k2bp_t* b, k2bp_t* c);
static void reck2bp_mul(k2bp_traversal_t* pos_a, const k2bp_t* a, k2bp_traversal_t* pos_b, const k2bp_t* b, k2bp_t* c);
static void reck2bp_scanmul(k2bp_traversal_t* pos_a, const k2bp_t* a, k2bp_traversal_t* pos_b, const k2bp_t* b, k2bp_t* c);

static uint8_t count(const uint64_t num);
static uint8_t count_p(const uint64_t num);
static uint64_t rank_p(const k2bp_traversal_t* pos_a, const k2bp_t* a);
static void build_rank_p(k2bp_traversal_t* pos_a, const k2bp_t* a);

static void k2bp_excdfs(k2bp_traversal_t* pos_a, const k2bp_t* a);
static void k2bp_scandfs(k2bp_traversal_t* pos_a, const k2bp_t* a);
static void k2bp_split(const k2bp_traversal_t* pos_a, const k2bp_t* a, k2bp_traversal_t* splits);
static void k2bp_splitinfo(const k2bp_traversal_t* pos_a, const k2bp_t* a, k2bp_traversal_t* splits);
static void k2bp_init_traversalinfo(const k2bp_traversal_t* pos_a, k2bp_traversal_t* copy_a);

// traverse nd copy subtree in a to c
static void k2bp_scandfs_copy(k2bp_traversal_t* pos_a, const k2bp_t* a, k2bp_t* c);
static void k2bp_excdfs_copy(k2bp_traversal_t* pos_a, const k2bp_t* a, k2bp_t* c);

static void reck2bp_decompress_subtrees(k2bp_traversal_t* pos_c, const k2bp_t* c, k2bp_t* a);

uint8_t k2bp_equal(const k2bp_t* a, const k2bp_t* b) {
  if(a->t.n != b->t.n) return 0;

  for(size_t i = 0; i < a->t.n; i++) {
    if(bv_i(&(a->t), i) != bv_i(&(b->t), i)) {
      return 0;
    }
  }

  if(a->n_l != b->n_l) return 0;
  for(size_t i = 0; i < a->n_l; i++) {
    if(k2bp_read_leaf(a, i) != k2bp_read_leaf(b, i)) {
      return 0;
    }
  }
  return 1;
}

// build k2 tree with BP representation
// from file `fname`, file has to be a text file with format
//  x1 y1
//  x2 y2
//  ...
//  xm ym
// where `m` is the amount of ones in the matrix
size_t k2bp_build_from_textfile(k2bp_t *a, const char* fname, size_t xsize) {
  assert(a != NULL && fname != NULL);
  FILE* f = fopen(fname, "rt");
  if(f == NULL) quit("k2bp_build_from_file: cannot open input file", __LINE__, __FILE__);

  size_t n; // number of entries
  size_t rmsize; // computed matrix size

  uint64_t *ia = create_ia(f,&n,&rmsize,xsize);
  assert(xsize==0 || rmsize==xsize);
  fclose(f);

  size_t size_tree = create_k2bp(ia, n, rmsize, a);
  free(ia);
  return size_tree;
}

// write a leaf in array `l`
// if there is not enough memory array is expanded by two
void k2bp_write_leaf(k2bp_t *a, uint8_t leaf){
  if(a->n_l >= a->maxn_l) {
    a->maxn_l *= 2;
    a->l = (uint8_t*) realloc(a->l, sizeof(uint8_t) * ((a->maxn_l + 1) / 2));
  }

  if(a->n_l % 2 == 0)
    a->l[a->n_l / 2] = leaf;
  else
    a->l[a->n_l / 2] = (leaf << 4) | (a->l[a->n_l / 2] & 15);
  a->n_l++;
}

// read leaf (4bits) from array `l`
uint8_t k2bp_read_leaf(const k2bp_t*a, const size_t pos) {
  assert(pos < a->n_l);
  if(pos % 2 == 0)
    return a->l[pos / 2] & 15;
  return (a->l[pos / 2] >> 4) & 15;
}

size_t k2bp_compute_height(const size_t rmsize) {
  assert(rmsize > 1);
  size_t msize = 2;
  // task for future me:
  //  add detection of overflow
  while(msize < rmsize) msize *= 2;
  return msize;
}

// free memory of a k2bp_t
void k2bp_free(k2bp_t* a) {
  a->msize = a->rmsize = a->m = a->maxn_l = a->n_l = a->n_p = 0;
  bv_free(&(a->t));
  if(a->l != NULL)
    free(a->l);

  if(a->exc_min_samples != NULL)
    free(a->exc_min_samples);
  if(a->exc_samples != NULL)
    free(a->exc_samples);
  if(a->leaves_samples != NULL)
    free(a->leaves_samples);

  if(a->subtreeinfo != NULL)
    free(a->subtreeinfo);
  if(a->leavesinfo != NULL)
    free(a->leavesinfo);

  if(a->pointers != NULL)
    free(a->pointers);

  a->l = NULL;
  a->exc_min_samples = NULL;
  a->exc_samples = NULL;
  a->leaves_samples = NULL;
  a->subtreeinfo = NULL;
  a->leavesinfo = NULL;
}

// return a pointer to an array with the coordinate of the non zeros
// the entries of the array the following format
//  x1 y1 x2 y2 ... xm ym
// n is the size of the array (by consequence n = 2 * m)
uint32_t* k2bp_nonzeros(const k2bp_t* a, size_t* n) {
  uint32_t* arr = (uint32_t*) malloc(sizeof(uint32_t) * (a->m * 2));
  *n = 0;
  k2bp_traversal_t pos_a = {a->msize, 0, 0, 0, 0, 0, 0, 0, 0};
  reck2bp_nonzeros(&pos_a, a, arr, n);
  return arr;
}

void k2bp_sum(const k2bp_t *a, const k2bp_t *b, k2bp_t *c) {
  assert(a != NULL && b != NULL && c != NULL);
  c->msize = a->msize;
  c->rmsize = a->rmsize;

  bv_init(&(c->t));
  c->maxn_l = 10;
  c->n_l = 0;
  c->l = (uint8_t*) malloc(sizeof(uint8_t) * c->maxn_l);
  c->m = 0;

  k2bp_traversal_t pos_a = K2BP_TRAVERSAL_INITIALIZER;
  k2bp_traversal_t pos_b = K2BP_TRAVERSAL_INITIALIZER;

  pos_a.msize = a->msize;
  pos_b.msize = b->msize;
  pos_a.excess = 1;
  pos_b.excess = 1;
  pos_a.size = a->t.n / 2;
  pos_b.size = b->t.n / 2;
  pos_a.leaves = a->n_l;
  pos_b.leaves = b->n_l;
  reck2bp_sum(&pos_a, a, &pos_b, b, c);
  return;
}

void k2bp_scansum(const k2bp_t *a, const k2bp_t *b, k2bp_t *c) {
  assert(a != NULL && b != NULL && c != NULL);
  c->msize = a->msize;
  c->rmsize = a->rmsize;

  bv_init(&(c->t));
  c->maxn_l = 10;
  c->n_l = 0;
  c->l = (uint8_t*) malloc(sizeof(uint8_t) * c->maxn_l);
  c->m = 0;

  k2bp_traversal_t pos_a = K2BP_TRAVERSAL_INITIALIZER;
  k2bp_traversal_t pos_b = K2BP_TRAVERSAL_INITIALIZER;

  pos_a.msize = a->msize;
  pos_b.msize = b->msize;
  pos_a.excess = 1;
  pos_b.excess = 1;
  pos_a.size = a->t.n / 2;
  pos_b.size = b->t.n / 2;
  pos_a.leaves = a->n_l;
  pos_b.leaves = b->n_l;
  reck2bp_scansum(&pos_a, a, &pos_b, b, c);
  return;
}

void k2bp_mul(const k2bp_t *a, const k2bp_t *b, k2bp_t *c) {
  assert(a != NULL && b != NULL && c != NULL);

  c->msize = a->msize;
  c->rmsize = a->rmsize;

  bv_init(&(c->t));
  c->maxn_l = 10;
  c->n_l = 0;
  c->l = (uint8_t*) malloc(sizeof(uint8_t) * c->maxn_l);
  c->m = 0;

  k2bp_traversal_t pos_a = K2BP_TRAVERSAL_INITIALIZER;
  k2bp_traversal_t pos_b = K2BP_TRAVERSAL_INITIALIZER;

  pos_a.msize = a->msize;
  pos_b.msize = b->msize;
  pos_a.excess = 1;
  pos_b.excess = 1;
  pos_a.size = a->t.n / 2;
  pos_b.size = b->t.n / 2;
  pos_a.leaves = a->n_l;
  pos_b.leaves = b->n_l;
  reck2bp_mul(&pos_a, a, &pos_b, b, c);
  return;
}

// save k2bp_t in files with prefix `fname`
// it creates possibly 4 files
//  `fname.i`  : stores just three integers that are the information of the matrix
//               `msize`, `rmsize` and `m`
//  `fname.t`  : stores the bit vector that represents the tree
//  `fname.l`  : stores the array with the leaves of the last level
//  `fname.exc`: stores the sampling information for a more efficient traversal
//               of smalls trees, if sampling is null the file is not created
//  `fname.r`  : stores the information of subtree information for a more
//               efficient traversal of big subtrees
void k2bp_save_to_file(const k2bp_t* a, const char* fname) {
  assert(a != NULL);
  assert(a->l != NULL && a->t.a != NULL);

  char info_name[1000];
  strcpy(info_name, fname);
  strcat(info_name, ".i");

  FILE* f = fopen(info_name, "w");
  if(f == NULL)
    quit("k2bp_save_to_file: file cannot be open", __LINE__, __FILE__);

  size_t w = fwrite(&(a->msize), sizeof(size_t), 1, f);
  if(w != 1)
    quit("k2bp_save_to_file: error writing in file", __LINE__, __FILE__);
  w = fwrite(&(a->rmsize), sizeof(size_t), 1, f);
  if(w != 1)
    quit("k2bp_save_to_file: error writing in file", __LINE__, __FILE__);
  w = fwrite(&(a->m), sizeof(size_t), 1, f);
  if(w != 1)
    quit("k2bp_save_to_file: error writing in file", __LINE__, __FILE__);

  fclose(f);

  char tree_name[1000];
  strcpy(tree_name, fname);
  strcat(tree_name, ".t");
  bv_save_to_file(&(a->t), tree_name);

  char l_name[1000];
  strcpy(l_name, fname);
  strcat(l_name, ".l");

  f = fopen(l_name, "w");
  if(f == NULL)
    quit("k2bp_save_to_file: file cannot be open", __LINE__, __FILE__);

  w = fwrite(&(a->maxn_l), sizeof(size_t), 1, f);
  if(w != 1)
    quit("k2bp_save_to_file: error writing in file", __LINE__, __FILE__);
  w = fwrite(&(a->n_l), sizeof(size_t), 1, f);
  if(w != 1)
    quit("k2bp_save_to_file: error writing in file", __LINE__, __FILE__);
  w = fwrite(a->l, sizeof(uint8_t), (a->n_l + 1)/ 2, f);
  if(w != (a->n_l + 1) / 2)
    quit("k2bp_save_to_file: error writing in file", __LINE__, __FILE__);

  fclose(f);

  if(a->exc_min_samples != NULL) {
    char exc_name[1000];
    strcpy(exc_name, fname);
    strcat(exc_name, ".exc");

    f = fopen(exc_name, "w");
    if(f == NULL)
      quit("k2bp_save_to_file: file cannot be open", __LINE__, __FILE__);

    w = fwrite(a->exc_min_samples, sizeof(uint16_t), SAMPLE_SIZE(a->t.n), f);
    if(w != SAMPLE_SIZE(a->t.n))
      quit("k2bp_save_to_file: error writing in file", __LINE__, __FILE__);
    w = fwrite(a->exc_samples, sizeof(uint16_t), SAMPLE_SIZE(a->t.n), f);
    if(w != SAMPLE_SIZE(a->t.n))
      quit("k2bp_save_to_file: error writing in file", __LINE__, __FILE__);
    w = fwrite(a->leaves_samples, sizeof(uint8_t), SAMPLE_SIZE(a->t.n), f);
    if(w != SAMPLE_SIZE(a->t.n))
      quit("k2bp_save_to_file: error writing in file", __LINE__, __FILE__);

    fclose(f);
  }

  if(a->subtreeinfo != NULL) {
    char sub_name[1000];
    strcpy(sub_name, fname);
    strcat(sub_name, ".s");

    f = fopen(sub_name, "w");
    if(f == NULL)
      quit("k2bp_save_to_file: file cannot be open", __LINE__, __FILE__);

    w = fwrite(&(a->n_info), sizeof(size_t), 1, f);
    if(w != 1)
      quit("k2bp_save_to_file: error writing in file", __LINE__, __FILE__);

    w = fwrite(&(a->threshold), sizeof(size_t), 1, f);
    if(w != 1)
      quit("k2bp_save_to_file: error writing in file", __LINE__, __FILE__);

    w = fwrite(a->subtreeinfo, sizeof(uint64_t), a->n_info, f);
    if(w != a->n_info)
      quit("k2bp_save_to_file: error writing in file", __LINE__, __FILE__);
    w = fwrite(a->leavesinfo, sizeof(uint32_t), a->n_info, f);
    if(w != a->n_info)
      quit("k2bp_save_to_file: error writing in file", __LINE__, __FILE__);

    fclose(f);
  }

  if(a->pointers != NULL) {
    char sub_name[1000];
    strcpy(sub_name, fname);
    strcat(sub_name, ".p");

    f = fopen(sub_name, "w");
    if(f == NULL)
      quit("k2bp_save_to_file: file cannot be open", __LINE__, __FILE__);

    w = fwrite(&(a->n_p), sizeof(size_t), 1, f);
    if(w != 1)
      quit("k2bp_save_to_file: error writing in file", __LINE__, __FILE__);
    w = fwrite(a->pointers, sizeof(size_t), a->n_p, f);
    if(w != a->n_p)
      quit("k2bp_save_to_file: error writing in file", __LINE__, __FILE__);
  }
}

// load k2bp_t from files with prefix `fname`
// it creates possibly 4 files
//  `fname.i`  : load just three integers that are the information of the matrix
//               `msize`, `rmsize` and `m`
//  `fname.t`  : load the bit vector that represents the tree
//  `fname.l`  : load the array with the leaves of the last level
//  `fname.exc`: load the sampling information for a more efficient traversal
//               of smalls trees, if sampling is null the file is not created
//  `fname.r`  : load the information of subtree information for a more
//               efficient traversal of big subtrees
// NOTE:
//   if `fname.(i|t|l)`
//    doesn't open, this function will throw an error
//   if `fname.(exc|r)`
//    doesn't open, doesn't matter means that the tree
//    doesn't has this extra information
void k2bp_load_from_file(k2bp_t* a, const char* fname) {
  assert(a != NULL);

  char info_name[1000];
  strcpy(info_name, fname);
  strcat(info_name, ".i");

  FILE* f = fopen(info_name, "r");
  if(f == NULL)
    quit("k2bp_load_from_file: file cannot be open", __LINE__, __FILE__);

  size_t w = fread(&(a->msize), sizeof(size_t), 1, f);
  if(w != 1)
    quit("k2bp_load_from_file: error reading from file", __LINE__, __FILE__);
  w = fread(&(a->rmsize), sizeof(size_t), 1, f);
  if(w != 1)
    quit("k2bp_load_from_file: error reading from file", __LINE__, __FILE__);
  w = fread(&(a->m), sizeof(size_t), 1, f);
  if(w != 1)
    quit("k2bp_load_from_file: error reading from file", __LINE__, __FILE__);

  fclose(f);

  char tree_name[1000];
  strcpy(tree_name, fname);
  strcat(tree_name, ".t");
  bv_load_from_file(&(a->t), tree_name);

  char l_name[1000];
  strcpy(l_name, fname);
  strcat(l_name, ".l");

  f = fopen(l_name, "r");
  if(f == NULL)
    quit("k2bp_load_from_file: file cannot be open", __LINE__, __FILE__);

  w = fread(&(a->maxn_l), sizeof(size_t), 1, f);
  if(w != 1)
    quit("k2bp_load_from_file: error reading from file", __LINE__, __FILE__);
  w = fread(&(a->n_l), sizeof(size_t), 1, f);
  if(w != 1)
    quit("k2bp_load_from_file: error reading from file", __LINE__, __FILE__);

  a->l = (uint8_t*) malloc(sizeof(uint8_t) * a->maxn_l / 2);

  w = fread(a->l, sizeof(uint8_t), (a->n_l + 1) / 2, f);
  if(w != (a->n_l + 1) / 2)
    quit("k2bp_load_from_file: error reading from file", __LINE__, __FILE__);

  fclose(f);

  char exc_name[1000];
  strcpy(exc_name, fname);
  strcat(exc_name, ".exc");

  f = fopen(exc_name, "r");
  if(f != NULL) {
    a->exc_min_samples = (uint16_t*) malloc(sizeof(uint16_t) * SAMPLE_SIZE(a->t.n));
    w = fread(a->exc_min_samples, sizeof(uint16_t), SAMPLE_SIZE(a->t.n), f);
    if(w != SAMPLE_SIZE(a->t.n))
      quit("k2bp_load_from_file: error reading from file", __LINE__, __FILE__);

    a->exc_samples = (uint16_t*) malloc(sizeof(uint16_t) * SAMPLE_SIZE(a->t.n));
    w = fread(a->exc_samples, sizeof(uint16_t), SAMPLE_SIZE(a->t.n), f);
    if(w != SAMPLE_SIZE(a->t.n))
      quit("k2bp_load_from_file: error reading from file", __LINE__, __FILE__);

    a->leaves_samples = (uint8_t*) malloc(sizeof(uint8_t) * SAMPLE_SIZE(a->t.n));
    w = fread(a->leaves_samples, sizeof(uint8_t), SAMPLE_SIZE(a->t.n), f);
    if(w != SAMPLE_SIZE(a->t.n))
      quit("k2bp_load_from_file: error reading from file", __LINE__, __FILE__);

    fclose(f);
  }

  char sub_name[100];
  strcpy(sub_name, fname);
  strcat(sub_name, ".s");

  f = fopen(sub_name, "r");
  if(f != NULL) {

    w = fread(&(a->n_info), sizeof(size_t), 1, f);
    if(w != 1)
      quit("k2bp_load_from_file: error reading from file", __LINE__, __FILE__);
    w = fread(&(a->threshold), sizeof(size_t), 1, f);
    if(w != 1)
      quit("k2bp_load_from_file: error reading from file", __LINE__, __FILE__);

    a->subtreeinfo = (uint64_t*) malloc(sizeof(uint64_t) * a->n_info);
    w = fread(a->subtreeinfo, sizeof(uint64_t), a->n_info, f);
    if(w != a->n_info)
      quit("k2bp_load_from_file: error reading from file", __LINE__, __FILE__);

    a->leavesinfo = (uint32_t*) malloc(sizeof(uint32_t) * a->n_info);
    w = fread(a->leavesinfo, sizeof(uint32_t), a->n_info, f);
    if(w != a->n_info)
      quit("k2bp_load_from_file: error reading from file", __LINE__, __FILE__);

    fclose(f);
  }

  char point_name[1000];
  strcpy(point_name, fname);
  strcat(point_name, ".p");

  f = fopen(point_name, "r");

  if(f != NULL) {
    w = fread(&(a->n_p), sizeof(size_t), 1, f);
    if(w != 1)
      quit("k2bp_save_to_file: error reading in file", __LINE__, __FILE__);
    a->pointers = (size_t*) malloc(sizeof(size_t) * a->n_p);
    w = fread(a->pointers, sizeof(size_t), a->n_p, f);
    if(w != a->n_p)
      quit("k2bp_save_to_file: error reading in file", __LINE__, __FILE__);
  }
}

void k2bp_dfs(k2bp_traversal_t* pos_a, const k2bp_t* a, size_t* nodes, size_t* leaves, size_t* nz, size_t* levels, size_t curr_level) {
  assert(pos_a->i_t < a->t.n && pos_a->i_l <= a->n_l);
  assert(bv_i(&(a->t), pos_a->i_t) == 1);
  curr_level++;
  if(pos_a->flag_p == 0)
    (*nodes)++;
  if(curr_level > *levels) *levels = curr_level;

  if(pos_a->i_t + 1 < a->t.n && bv_i(&(a->t), pos_a->i_t + 1) == 0) {
    pos_a->i_t += 2;
    return;
  }
  
  if(pos_a->msize == _K_) {
    uint8_t leaf = k2bp_read_leaf(a, pos_a->i_l);
    pos_a->i_l++;
    pos_a->i_t += 4;
    if(pos_a->flag_p == 0)
      (*leaves)++;

    if(leaf & 1) {
      (*nz)++;
    }
    if(leaf & 2) {
      (*nz)++;
    }
    if(leaf & 4) {
      (*nz)++;
    }
    if(leaf & 8) {
      (*nz)++;
    }
    return;
  }

  if(bv_get_int(&(a->t), pos_a->i_t, 6) == LEAF_P) {
    size_t det = a->pointers[pos_a->i_p];
    k2bp_traversal_t new_pos = K2BP_TRAVERSAL_INITIALIZER;
    k2bp_init_traversalinfo(pos_a, &(new_pos));
    new_pos.i_t = det;
    new_pos.i_l = pos_a->i_l;
    new_pos.i_p = rank_p(&new_pos, a);
    new_pos.flag_p = 1;
    k2bp_dfs(&new_pos, a, nodes, leaves, nz, levels, curr_level);
    if(pos_a->flag_p == 0)
      *(nodes) += 2;
    pos_a->i_l = new_pos.i_l;
    pos_a->i_t += 6;
    pos_a->i_p++;
    return;
  }

  k2bp_traversal_t pos_aux = {pos_a->msize / 2, pos_a->x, pos_a->y, pos_a->i_t + 1, pos_a->i_l, 0, 0, 0, 0, 0, 0, NULL};
  pos_aux.i_p = pos_a->i_p;
  pos_aux.rank_pointers = pos_a->rank_pointers;
  pos_aux.flag_p = pos_a->flag_p;
  k2bp_dfs(&pos_aux, a, nodes, leaves, nz, levels, curr_level);

  pos_aux.x = pos_a->x;
  pos_aux.y = pos_a->y + pos_a->msize / 2;
  pos_aux.flag_p = pos_a->flag_p;
  k2bp_dfs(&pos_aux, a, nodes, leaves, nz, levels, curr_level);

  pos_aux.x = pos_a->x + pos_a->msize / 2;
  pos_aux.y = pos_a->y;
  pos_aux.flag_p = pos_a->flag_p;
  k2bp_dfs(&pos_aux, a, nodes, leaves, nz, levels, curr_level);

  pos_aux.x = pos_a->x + pos_a->msize / 2;
  pos_aux.y = pos_a->y + pos_a->msize / 2;
  pos_aux.flag_p = pos_a->flag_p;
  k2bp_dfs(&pos_aux, a, nodes, leaves, nz, levels, curr_level);
  pos_a->i_t = pos_aux.i_t + 1;
  pos_a->i_l = pos_aux.i_l;
  pos_a->i_p = pos_aux.i_p;
  
}

size_t k2bp_stats(const k2bp_t* a, size_t* nodes, size_t* leaves, size_t* nz) {
  size_t levels = 0;
  *nodes = *nz = *leaves = 0;
  k2bp_traversal_t pos_a = K2BP_TRAVERSAL_INITIALIZER;
  pos_a.msize = a->msize;
  if(a->pointers != NULL)
    build_rank_p(&pos_a, a);
  k2bp_dfs(&pos_a, a, nodes, leaves, nz, &levels, 0);
  if(pos_a.rank_pointers != NULL)
    free(pos_a.rank_pointers);
  return levels;
}

size_t k2bp_show_stats(const k2bp_t *a, const char *fname, FILE *f) {
  fprintf(f, "file: %s\n", fname);
  fprintf(f, "matrix size: %zu, leaf size: %d, k2 internal size: %zu\n", a->rmsize, _K_, a->msize);

  size_t nodes, leaves, nz;
  nodes = leaves = nz = 0;
  size_t levels = k2bp_stats(a, &nodes, &leaves, &nz);
  assert(nz == a->m);
  assert(((nodes + leaves) * 2) == a->t.n);

  fprintf(f, " nonzeros: %zu, nonzeros x row: %.3lf\n", nz, (double) nz / a->rmsize);
  fprintf(f, " levels: %zu, nodes: %zu, leaves: %zu\n", levels, nodes, leaves);
  size_t bp_bytes = sizeof(bv_t) + (a->t.n + 64 - 1) / 64 * sizeof(uint64_t);
  fprintf(f, " size by parts\n");
  fprintf(f, "  bp  size: %zu bytes, %zu bits, %.3lf bits x nonzero\n", bp_bytes, bp_bytes * CHAR_BIT, (double) bp_bytes * CHAR_BIT / nz);
  size_t l_bytes = sizeof(size_t) * 2 + sizeof(uint8_t) * (a->n_l + 1) / 2;
  fprintf(f, "  l   size: %zu bytes, %zu bits, %.3lf bits x nonzero\n", l_bytes, l_bytes * CHAR_BIT, (double) l_bytes * CHAR_BIT / nz);
  size_t exc_bytes = 0;
  size_t mic_bytes = 0;
  if(a->exc_min_samples != NULL) {
    exc_bytes = SAMPLE_SIZE(a->t.n) * (sizeof(uint16_t) * 2 + sizeof(uint8_t));
    mic_bytes = 65536 * sizeof(uint8_t) * 2;
  }
  fprintf(f, "  exc size: %zu bytes, %zu bits, %.3lf bits x nonzero\n", exc_bytes, exc_bytes * CHAR_BIT, (double) exc_bytes * CHAR_BIT / nz);
  fprintf(f, "  mic size: %zu bytes, %zu bits, %.3lf bits x nonzero\n", mic_bytes, mic_bytes * CHAR_BIT, (double) mic_bytes * CHAR_BIT / nz);
  size_t sub_bytes = 0;
  if(a->subtreeinfo != NULL) {
    sub_bytes = a->n_info * (sizeof(uint64_t) + sizeof(uint32_t)) + sizeof(size_t) * 2;
  }
  fprintf(f, "  sub size: %zu bytes, %zu bits, %.3lf bits x nonzero\n", sub_bytes, sub_bytes * CHAR_BIT, (double) sub_bytes * CHAR_BIT / nz);
  size_t pointer_bytes = 0;
  if(a->pointers != NULL) {
    pointer_bytes = a->n_p * sizeof(size_t);
  }
  fprintf(f, "  ptr size: %zu bytes, %zu bits, %.3lf bits x nonzero\n", pointer_bytes, pointer_bytes * CHAR_BIT, (double) pointer_bytes * CHAR_BIT / nz);
  size_t total_bytes = bp_bytes + l_bytes + exc_bytes + sub_bytes + mic_bytes + sizeof(size_t) * 3 + pointer_bytes;
  fprintf(f, " total size: %zu bytes, %zu bits, %.3lf bits x nonzero\n", total_bytes, total_bytes * CHAR_BIT, (double) total_bytes * CHAR_BIT / nz);

  printf("%zu\n", a->n_p);
  size_t max = 0;
  for(size_t i = 0; i < a->n_p; i++) {
    if(ceil_log2(a->pointers[i]) > max) max = ceil_log2(a->pointers[i]);
  }
  printf("%zu\n", max);
  printf("%.3lf\n", (double) max * a->n_p / nz);
  printf("%.3lf\n", (double) (total_bytes * CHAR_BIT - pointer_bytes * CHAR_BIT + max * a->n_p)/nz);

  iv_t z;
  iv_init(&z, a->n_p, max);
  for(size_t i = 0; i < a->n_p; i++) {
    iv_set(&z, i, a->pointers[i]);
  }
  for(size_t i = 0; i < a->n_p; i++) {
    assert(iv_get(&z, i) == a->pointers[i]);
  }
  return total_bytes;
}

// build excess sampling
// this helps for small trees across `a`
void k2bp_build_exc_sampling(k2bp_t* a) {
  assert(a != NULL);
  assert(a->t.a != NULL && a->l != NULL);
  if(a->exc_samples != NULL) { // clean what is there
    free(a->exc_samples);
    free(a->exc_min_samples);
    free(a->leaves_samples);
    a->exc_min_samples = a->exc_samples = NULL;
    a->leaves_samples = NULL;
  }

  a->exc_samples = (uint16_t*) malloc(sizeof(uint16_t) * SAMPLE_SIZE(a->t.n));
  a->exc_min_samples = (uint16_t*) malloc(sizeof(uint16_t) * SAMPLE_SIZE(a->t.n));
  a->leaves_samples = (uint8_t*) malloc(sizeof(uint8_t) * SAMPLE_SIZE(a->t.n));
  uint16_t excess = 1;
  uint64_t min_excess = 1;
  uint8_t leaves = 0;
  size_t block = 0;
  for(size_t i = 1; i < a->t.n; i++) {
    if(i % BLOCK_SIZE == 0) {
      a->exc_samples[block] = excess;
      a->exc_min_samples[block] = min_excess;
      a->leaves_samples[block] = leaves;
      block++;
      min_excess = LLONG_MAX;
      leaves = 0;
    }
    if(i + 3 < a->t.n) {
      if(bv_get_int(&(a->t), i, 4) == LEAF_1) leaves++;
    }
    excess += ((bv_i(&(a->t), i)) ? 1 : -1);
    if(excess < min_excess)
      min_excess = excess;
  }
  a->exc_samples[block] = excess;
  a->exc_min_samples[block] = min_excess;
  a->leaves_samples[block] = leaves;
}

void k2bp_addsubtree_info(k2bp_t* a, size_t threshold) {
  assert(a != NULL);
  assert(a->t.a != NULL && a->l != NULL);
  a->threshold = threshold;
  if(a->subtreeinfo != NULL) {
    free(a->subtreeinfo);
    free(a->leavesinfo);
    a->subtreeinfo = NULL; a->leavesinfo = NULL;
    a->n_info = 0;
  }

  k2bp_traversal_t pos_a = {a->msize, 0, 0, 0, 0, 0, 0, 0, 0};

  uint32_t leaves = 0;
  vu64_t subinfo, leavesinfo;
  vu64_init(&subinfo);
  vu64_init(&leavesinfo);

  reck2bp_addsubtree_info(&pos_a, a, &leaves, &subinfo, &leavesinfo);

  assert(subinfo.n == leavesinfo.n);

  a->subtreeinfo = (uint64_t*) malloc(sizeof(uint64_t) * subinfo.n);
  a->leavesinfo = (uint32_t*) malloc(sizeof(uint32_t) * leavesinfo.n);
  a->n_info = subinfo.n;

  for(size_t i = 0; i < a->n_info; i++) {
    a->subtreeinfo[i] = subinfo.v[i];
    a->leavesinfo[i] = leavesinfo.v[i];
  }

  vu64_free(&subinfo);
  vu64_free(&leavesinfo);
}

size_t k2bp_checksubtree_info(const k2bp_t* a) {
  assert(a != NULL);
  assert(a->t.a != NULL && a->l != NULL && a->subtreeinfo != NULL && a->leavesinfo != NULL);

  k2bp_traversal_t pos_a = {0, 0, 0, 0, 0, 0, a->t.n / 2, a->n_l, 0};
  reck2bp_checksubtree_info(&pos_a, a);
  return 1;
}

void k2bp_compress_subtrees(const k2bp_t* a, k2bp_t* c, size_t limit) {
  assert(a != NULL && c != NULL);
  assert(limit > 0);

  c->msize = a->msize;
  c->rmsize = a->rmsize;
  c->m = a->m;

  uint8_t* text = (uint8_t*) malloc(sizeof(uint8_t) * a->t.n);
  int64_t* csa = (int64_t*) malloc(sizeof(int64_t) * a->t.n);
  int64_t* plcp = (int64_t*) malloc(sizeof(int64_t) * a->t.n);
  int64_t* lcp = (int64_t*) malloc(sizeof(int64_t) * a->t.n);

  for(size_t i = 0; i < a->t.n; i++) {
    text[i] = bv_i(&(a->t), i);
    assert(text[i] == 1 || text[i] == 0);
  }

  if(libsais64(text, csa, a->t.n, 0, NULL) != 0)
    quit("k2bp_compress_subtrees: error building csa", __LINE__, __FILE__);
  if(libsais64_plcp(text, csa, plcp, a->t.n) != 0)
    quit("k2bp_compress_subtrees: error building plcp", __LINE__, __FILE__);
  if(libsais64_lcp(plcp, csa, lcp, a->t.n) != 0)
    quit("k2bp_compress_subtrees: error building plcp", __LINE__, __FILE__);


  dsu groups;
  dsu_init(&groups, a->t.n);

  for(size_t i = 1; i < a->t.n; i++) { // This can be optimized to t.n/2 I think
                                       // test it later
    if(bv_i(&(a->t), csa[i]) == 0) continue; // not real tree
    if(bv_i(&(a->t), csa[i - 1]) == 0) continue; // cannot compared against no real tree

    k2bp_traversal_t pos_a_i = K2BP_TRAVERSAL_INITIALIZER;
    k2bp_traversal_t pos_a_i_1 = K2BP_TRAVERSAL_INITIALIZER;
    pos_a_i.i_t = csa[i]; pos_a_i.excess = 1;
    pos_a_i_1.i_t = csa[i - 1]; pos_a_i_1.excess = 1;
    k2bp_scandfs(&pos_a_i, a);
    k2bp_scandfs(&pos_a_i_1, a);

    if(lcp[i] >= pos_a_i.size * 2) { // found identical subtrees
      assert(pos_a_i.size == pos_a_i_1.size);
      dsu_union_set(&groups, csa[i], csa[i - 1]);
    }
  }

  free(text);
  free(plcp);
  free(csa);
  free(lcp);

  uint64_t* prefix_help = (uint64_t*) malloc(sizeof(uint64_t) * a->t.n);
  for(size_t i = 0; i < a->t.n; i++)
    prefix_help[i] = 0;
  vu64_t pointers;
  vu64_init(&pointers);
  bv_init(&(c->t));

  for(size_t i = 0; i < a->t.n; i++) {
    if(bv_i(&(a->t), i)) {
      bv_pb(&(c->t), 1);
      size_t left = dsu_find_set(&groups, i);
      k2bp_traversal_t pos_a_i = K2BP_TRAVERSAL_INITIALIZER;
      pos_a_i.i_t = i;
      pos_a_i.excess = 1;
      k2bp_scandfs(&pos_a_i, a);
      
      if(left != i && limit + 6 <= pos_a_i.size * 2) {
        vu64_grow(&pointers, 1);
        pointers.v[pointers.n - 1] = left - prefix_help[left - 1];
        bv_pb(&(c->t), 1);
        bv_pb(&(c->t), 0);
        bv_pb(&(c->t), 1);
        bv_pb(&(c->t), 0);
        bv_pb(&(c->t), 0);

        for(size_t j = i; j < i + pos_a_i.size * 2; j++) {
          prefix_help[j] = prefix_help[j - 1] + 1;
        }
        prefix_help[i + pos_a_i.size * 2 - 1] -= 6;
        i = i + pos_a_i.size * 2 - 1;
      } else {
        if(i > 0) prefix_help[i] = prefix_help[i - 1];
      }
    } else {
      bv_pb(&(c->t), 0);
      prefix_help[i] = prefix_help[i - 1];
    }
  }
  bv_shrink(&(c->t));
  
  c->n_p = pointers.n;
  c->pointers = (size_t*) malloc(sizeof(size_t) * pointers.n);
  for(size_t i = 0; i < pointers.n; i++) {
    c->pointers[i] = pointers.v[i];
  }

  vu64_free(&pointers);

  free(prefix_help);
  c->n_l = 0;
  c->maxn_l = a->maxn_l;
  c->l = (uint8_t*) malloc(sizeof(uint8_t) * (a->maxn_l + 1) / 2);
  for(size_t i = 0; i < a->n_l; i++) {
    k2bp_write_leaf(c, k2bp_read_leaf(a, i));
  }

  if(a->subtreeinfo != NULL) {
    c->threshold = a->threshold;
    c->n_info = a->n_info;
    c->subtreeinfo = (uint64_t*) malloc(sizeof(uint64_t) * a->n_info);
    c->leavesinfo = (uint32_t*) malloc(sizeof(uint32_t) * a->n_info);
    for(size_t i = 0; i < a->n_info; i++) {
      c->subtreeinfo[i] = a->subtreeinfo[i];
      c->leavesinfo[i] = a->leavesinfo[i];
    }
  }

  if(a->exc_min_samples != NULL) {
    k2bp_build_exc_sampling(c);
  }
}

void k2bp_decompress_subtrees(const k2bp_t *c, k2bp_t *a) {
  a->msize = c->msize;
  a->rmsize = c->rmsize;
  a->m = c->m;

  a->maxn_l = c->maxn_l;
  a->n_l = 0;
  a->l = (uint8_t*) malloc(sizeof(uint8_t) * (a->maxn_l + 1) / 2);

  if(c->subtreeinfo != NULL) {
    a->n_info = c->n_info;
    a->threshold = c->threshold;
    a->subtreeinfo = (uint64_t*) malloc(sizeof(uint64_t) * c->n_info);
    a->leavesinfo = (uint32_t*) malloc(sizeof(uint32_t) * c->n_info);
    for(size_t i = 0; i < c->n_info; i++) {
      a->subtreeinfo[i] = c->subtreeinfo[i];
      a->leavesinfo[i] = c->leavesinfo[i];
    }
  }

  bv_init(&(a->t));
  k2bp_traversal_t pos_c = K2BP_TRAVERSAL_INITIALIZER;
  pos_c.msize = c->msize;

  build_rank_p(&pos_c, c);

  reck2bp_decompress_subtrees(&pos_c, c, a);
  free(pos_c.rank_pointers);
  assert(a->n_l == c->n_l);
}

// ----------------------------------------------------------

// auxiliary functions
static void reck2bp_decompress_subtrees(k2bp_traversal_t* pos_c, const k2bp_t* c, k2bp_t* a) {
//  printf("%zu %zu %zu %zu\n", pos_c->msize, pos_c->i_t, pos_c->i_l, pos_c->i_p);
//  for(size_t i = pos_c->i_t; i < pos_c->i_t + 20; i++) {
//    if(bv_i(&(c->t), i)) printf("(");
//    else printf(")");
//  }
//  printf("\n");
  assert(bv_i(&(c->t), pos_c->i_t) == 1);

  if(bv_get_int(&(c->t), pos_c->i_t, 2) == LEAF_0) {
    bv_pb(&(a->t), 1);
    bv_pb(&(a->t), 0);

    pos_c->i_t += 2;
    return;
  }

  if(pos_c->msize == _K_) {
    bv_pb(&(a->t), 1);
    bv_pb(&(a->t), 1);
    bv_pb(&(a->t), 0);
    bv_pb(&(a->t), 0);
    k2bp_write_leaf(a, k2bp_read_leaf(c, pos_c->i_l));

    pos_c->i_t += 4;
    pos_c->i_l++;
    return;
  }

  if(bv_get_int(&(c->t), pos_c->i_t, 6) == LEAF_P) {
    size_t det = c->pointers[pos_c->i_p];
    k2bp_traversal_t new_pos = K2BP_TRAVERSAL_INITIALIZER;
    k2bp_init_traversalinfo(pos_c, &(new_pos));
    new_pos.i_t = det;
    new_pos.i_l = pos_c->i_l;
    new_pos.i_p = rank_p(&new_pos, c);
    new_pos.flag_p = 1;
    reck2bp_decompress_subtrees(&new_pos, c, a);
    pos_c->i_l = new_pos.i_l;
    pos_c->i_t += 6;
    pos_c->i_p++;
    return;
  }

  bv_pb(&(a->t), 1);
  size_t curr_msize = pos_c->msize;
  pos_c->i_t++;
  pos_c->msize = curr_msize / 2;
  reck2bp_decompress_subtrees(pos_c, c, a);
  pos_c->msize = curr_msize / 2;
  reck2bp_decompress_subtrees(pos_c, c, a);
  pos_c->msize = curr_msize / 2;
  reck2bp_decompress_subtrees(pos_c, c, a);
  pos_c->msize = curr_msize / 2;
  reck2bp_decompress_subtrees(pos_c, c, a);
  bv_pb(&(a->t), 0);
  pos_c->i_t++;
  pos_c->msize = curr_msize;
}

static uint8_t count(const uint64_t num) {
  uint64_t x = num;
  uint64_t y = num >> 1;
  uint64_t hi = x & y;
  uint64_t lo = ~ (x | y);

  uint64_t bits = (hi & (lo >> 2)) & 2305843009213693951;
  return __builtin_popcountll(bits);
}

static uint8_t count_p(const uint64_t num) {
  uint64_t bits = ((num & (num >> 1) & (~(num >> 2)) & (num >> 3) & (~(num >> 4)) & (~(num >> 5)))) & (((uint64_t)-1) >> 5);
  return __builtin_popcountll(bits);
}

static uint64_t rank_p(const k2bp_traversal_t* pos_a, const k2bp_t* a) {
  const uint64_t block = pos_a->i_t / BLOCK_SIZE;
  uint64_t ret = pos_a->rank_pointers[block];
  size_t curr_i = BLOCK_SIZE * (pos_a->i_t / BLOCK_SIZE);
  for(; curr_i + 64 < pos_a->i_t; curr_i += 64) {
    uint64_t bits = bv_get_int(&(a->t), curr_i, 64);
    ret += count_p(bits);
    if(curr_i + 64 < pos_a->i_t) {
      uint64_t len = (curr_i + 69 < a->t.n ? 10 :  a->t.n - (curr_i + 59));
      ret += count_p(bv_get_int(&(a->t), curr_i + 59, len) | (((uint64_t) -1) << len));
    }
  }

  for(; curr_i < pos_a->i_t; curr_i++) {
    if(curr_i + 6 > pos_a->i_t) break;
    if(bv_get_int(&(a->t), curr_i, 6) == LEAF_P) ret++;
  }
  return ret;
}

static void build_rank_p(k2bp_traversal_t* pos_a, const k2bp_t* a) {
  if(pos_a->rank_pointers != NULL) return;
  
  pos_a->rank_pointers = (uint64_t*) malloc(sizeof(uint64_t) * SAMPLE_SIZE(a->t.n));
  size_t block = 0;
  uint64_t prefix_sum = 0;
  size_t i;
  for(i = 0; i + 64 < a->t.n; i += 64) {
    if(i % BLOCK_SIZE == 0) {
      pos_a->rank_pointers[block] = prefix_sum;
      block++;
    }

    uint64_t bits = bv_get_int(&(a->t), i, 64);
    prefix_sum += count_p(bits);
    prefix_sum += count_p(bv_get_int(&(a->t), i + 59, 10) | (((uint64_t) -1) << 10));
  }

  for(; i + 6 < a->t.n; i++) {
    if(i % BLOCK_SIZE == 0) {
      pos_a->rank_pointers[block] = prefix_sum;
      block++;
    }
    prefix_sum += bv_get_int(&(a->t), i, 6) == LEAF_P;
  }

  pos_a->rank_pointers[block] = prefix_sum;
  block++;

}

static void k2bp_scandfs_copy(k2bp_traversal_t* pos_a, const k2bp_t* a, k2bp_t* c) {
  assert(bv_i(&(a->t), pos_a->i_t) == 1);
  if(bv_get_int(&(a->t), pos_a->i_t, 4) == LEAF_1) {
    bv_pb(&(c->t), 1);
    bv_pb(&(c->t), 1);
    bv_pb(&(c->t), 0);
    bv_pb(&(c->t), 0);
    pos_a->i_t += 4;
    pos_a->size = 2;
    pos_a->leaves = 1;
    k2bp_write_leaf(c, k2bp_read_leaf(a, pos_a->i_l));
    c->m += __builtin_popcount(k2bp_read_leaf(a, pos_a->i_l));
    pos_a->i_l++;
    return;
  }
  if(bv_get_int(&(a->t), pos_a->i_t, 2) == LEAF_0) {
    bv_pb(&(c->t), 1);
    bv_pb(&(c->t), 0);
    pos_a->i_t += 2;
    pos_a->size = 1;
    pos_a->leaves = 0;
    return;
  }

  pos_a->size = 0;
  pos_a->leaves = 0;
  int64_t obj_excess = pos_a->excess;
  uint64_t curr_pos = pos_a->i_t;
  bv_pb(&(c->t), 1);
  pos_a->i_t++;
  for(; pos_a->i_t + 16 <= a->t.n; pos_a->i_t += 16) {
    uint64_t bits = bv_get_int(&(a->t), pos_a->i_t, 16);
    // found microblock
    if(obj_excess >= pos_a->excess + exc_min_micro[bits] + 1) {
      for(uint8_t i = 0; i < 16; i++) {
        pos_a->i_t++;
        if(bits & (1ULL << i)) {
          pos_a->excess++;
          bv_pb(&(c->t), 1);
        } else {
          pos_a->excess--;
          bv_pb(&(c->t), 0);
        }
        if(obj_excess == pos_a->excess + 1) {
          pos_a->size = (pos_a->i_t - curr_pos) / 2;
          pos_a->leaves += count(bits | (-1ULL << (i + 1)));
          for(size_t i = 0; i < pos_a->leaves; i++) {
            c->m += __builtin_popcount(k2bp_read_leaf(a, pos_a->i_l + i));
            k2bp_write_leaf(c, k2bp_read_leaf(a, pos_a->i_l + i));
          }
          pos_a->i_l += pos_a->leaves;
          return;
        }
      }
    }
    bv_append_u16(&(c->t), bits);
    assert(bv_get_int(&(c->t), c->t.n - 16, 16) == bits);
    pos_a->excess += exc_micro[bits];
    if(pos_a->i_t + 19 <= a->t.n) {
      pos_a->leaves += count(bv_get_int(&(a->t), pos_a->i_t, 19) |
          (-1ULL << 19));
    } else {
      pos_a->leaves += count(bv_get_int(&(a->t), pos_a->i_t, (a->t.n - pos_a->i_t)) |
          (-1ULL << (a->t.n - pos_a->i_t)));
    }
  }
  for(;;) {
    if(pos_a->i_t + 4 <= a->t.n && bv_get_int(&(a->t), pos_a->i_t, 4) == LEAF_1) {
      pos_a->leaves++;
    }
    if(bv_i(&(a->t), pos_a->i_t)) {
      bv_pb(&(c->t), 1);
      pos_a->excess++;
    } else {
      bv_pb(&(c->t), 0);
      pos_a->excess--;
    }
    pos_a->i_t++;
    if(obj_excess == pos_a->excess + 1) {
      pos_a->size = (pos_a->i_t - curr_pos) / 2;
      for(size_t i = 0; i < pos_a->leaves; i++) {
        c->m += __builtin_popcount(k2bp_read_leaf(a, pos_a->i_l + i));
        k2bp_write_leaf(c, k2bp_read_leaf(a, pos_a->i_l + i));
      }
      pos_a->i_l += pos_a->leaves;
      return;
    }
  }
}

static void k2bp_scandfs(k2bp_traversal_t* pos_a, const k2bp_t* a) {
  assert(bv_i(&(a->t), pos_a->i_t) == 1);
  if(bv_get_int(&(a->t), pos_a->i_t, 4) == LEAF_1) {
    pos_a->i_t += 4;
    pos_a->size = 2;
    pos_a->leaves = 1;
    return;
  }
  if(bv_get_int(&(a->t), pos_a->i_t, 2) == LEAF_0) {
    pos_a->i_t += 2;
    pos_a->size = 1;
    pos_a->leaves = 0;
    return;
  }

  pos_a->size = 0;
  pos_a->leaves = 0;
  int64_t obj_excess = pos_a->excess;
  uint64_t curr_pos = pos_a->i_t;
  pos_a->i_t++;
  for(; pos_a->i_t + 16 <= a->t.n; pos_a->i_t += 16) {
    uint64_t bits = bv_get_int(&(a->t), pos_a->i_t, 16);
    // found microblock
    if(obj_excess >= pos_a->excess + exc_min_micro[bits] + 1) {
      for(uint8_t i = 0; i < 16; i++) {
        pos_a->i_t++;
        if(bits & (1ULL << i)) pos_a->excess++;
        else pos_a->excess--;
        if(obj_excess == pos_a->excess + 1) {
          pos_a->size = (pos_a->i_t - curr_pos) / 2;
          pos_a->leaves += count(bits | (-1ULL << (i + 1)));
          return;
        }
      }
    }
    pos_a->excess += exc_micro[bits];
    if(pos_a->i_t + 19 <= a->t.n) {
      pos_a->leaves += count(bv_get_int(&(a->t), pos_a->i_t, 19) |
          (-1ULL << 19));
    } else {
      pos_a->leaves += count(bv_get_int(&(a->t), pos_a->i_t, (a->t.n - pos_a->i_t)) |
          (-1ULL << (a->t.n - pos_a->i_t)));
    }
  }
  for(;;) {
    if(pos_a->i_t + 4 <= a->t.n && bv_get_int(&(a->t), pos_a->i_t, 4) == LEAF_1) {
      pos_a->leaves++;
    }
    if(bv_i(&(a->t), pos_a->i_t)) pos_a->excess++;
    else pos_a->excess--;
    pos_a->i_t++;
    if(obj_excess == pos_a->excess + 1) {
      pos_a->size = (pos_a->i_t - curr_pos) / 2;
      return;
    }
  }
}

static void k2bp_excdfs_copy(k2bp_traversal_t* pos_a, const k2bp_t* a, k2bp_t* c) {
  assert(bv_i(&(a->t), pos_a->i_t) == 1);
  if(bv_get_int(&(a->t), pos_a->i_t, 4) == LEAF_1) {
    bv_pb(&(c->t), 1);
    bv_pb(&(c->t), 1);
    bv_pb(&(c->t), 0);
    bv_pb(&(c->t), 0);
    pos_a->i_t += 4;
    pos_a->size = 2;
    pos_a->leaves = 1;
    k2bp_write_leaf(c, k2bp_read_leaf(a, pos_a->i_l));
    c->m += __builtin_popcount(k2bp_read_leaf(a, pos_a->i_l));
    pos_a->i_l++;
    return;
  }
  if(bv_get_int(&(a->t), pos_a->i_t, 2) == LEAF_0) {
    bv_pb(&(c->t), 1);
    bv_pb(&(c->t), 0);
    pos_a->i_t += 2;
    pos_a->size = 1;
    pos_a->leaves = 0;
    return;
  }

  pos_a->leaves = 0;
  pos_a->size = 0;
  size_t curr_pos = pos_a->i_t;
  int64_t obj_excess = pos_a->excess;
  bv_pb(&(c->t), 1);
  pos_a->i_t++;
  size_t obj_pos = BLOCK_SIZE * ((pos_a->i_t + BLOCK_SIZE - 1) / BLOCK_SIZE);
  // micro tables
  for(; pos_a->i_t + 16 < obj_pos; pos_a->i_t += 16) {
    uint64_t bits = bv_get_int(&(a->t), pos_a->i_t, 16);
    // found micro block
    if(obj_excess >= pos_a->excess + exc_min_micro[bits] + 1) {
      for(uint8_t i = 0; i < 16; i++) {
        pos_a->i_t++;
        if(bits & (1ULL << i)) {
          pos_a->excess++;
          bv_pb(&(c->t), 1);
        } else {
          pos_a->excess--;
          bv_pb(&(c->t), 0);
        }

        if(obj_excess == pos_a->excess + 1) {
          pos_a->size = (pos_a->i_t - curr_pos) / 2;
          pos_a->leaves += count(bits | (-1ULL << (i + 1)));
          for(size_t i = 0; i < pos_a->leaves; i++) {
            k2bp_write_leaf(c, k2bp_read_leaf(a, pos_a->i_l + i));
            c->m += __builtin_popcount(k2bp_read_leaf(a, pos_a->i_l + i));
          }
          pos_a->i_l += pos_a->leaves;
          return;
        }
      }
    }
    bv_append_u16(&(c->t), bits);
    pos_a->excess += exc_micro[bits];
    pos_a->leaves += count(bv_get_int(&(a->t), pos_a->i_t, 19) |
                           (-1ULL << 19));
  }
  if(pos_a->i_t % BLOCK_SIZE != 0) {
    size_t extra = BLOCK_SIZE - pos_a->i_t % BLOCK_SIZE;
    uint64_t bits = bv_get_int(&(a->t), pos_a->i_t, extra);
    size_t save_pos = pos_a->i_t;
    for(uint8_t i = 0; i < extra; i++) {
      pos_a->i_t++;
      if(bits & (1ULL << i)) {
        pos_a->excess++;
        bv_pb(&(c->t), 1);
      } else {
        pos_a->excess--;
        bv_pb(&(c->t), 0);
      }
      if(obj_excess == pos_a->excess + 1) {
        pos_a->size = (pos_a->i_t - curr_pos) / 2;
        pos_a->leaves += count(bits | (-1ULL << (i + 1)));
        for(size_t i = 0; i < pos_a->leaves; i++) {
          k2bp_write_leaf(c, k2bp_read_leaf(a, pos_a->i_l + i));
          c->m += __builtin_popcount(k2bp_read_leaf(a, pos_a->i_l + i));
        }
        pos_a->i_l += pos_a->leaves;
        return;
      }
    }
    pos_a->leaves += count(bv_get_int(&(a->t), save_pos, extra + 3) | (-1ULL << (extra + 3)));
  }

  assert(pos_a->i_t % BLOCK_SIZE == 0);

  // jumps between blocks
  size_t block = pos_a->i_t / BLOCK_SIZE;
  for(;; pos_a->i_t += BLOCK_SIZE) {
    // found block
    if(obj_excess > a->exc_min_samples[block]) {
      for(; pos_a->i_t + 16 <= a->t.n; pos_a->i_t += 16) {
        uint64_t bits = bv_get_int(&(a->t), pos_a->i_t, 16);
        // found microblock
        if(obj_excess >= pos_a->excess + exc_min_micro[bits] + 1) {
          for(uint8_t i = 0; i < 16; i++) {
            pos_a->i_t++;
            if(bits & (1ULL << i)) {
              pos_a->excess++;
              bv_pb(&(c->t), 1);
            } else {
              pos_a->excess--;
              bv_pb(&(c->t), 0);
            }
            if(obj_excess == pos_a->excess + 1) {
              pos_a->size = (pos_a->i_t - curr_pos) / 2;
              pos_a->leaves += count(bits | (-1ULL << (i + 1)));
              for(size_t i = 0; i < pos_a->leaves; i++) {
                k2bp_write_leaf(c, k2bp_read_leaf(a, pos_a->i_l + i));
                c->m += __builtin_popcount(k2bp_read_leaf(a, pos_a->i_l + i));
              }
              pos_a->i_l += pos_a->leaves;
              return;
            }
          }
        }
        bv_append_u16(&(c->t), bits);
        pos_a->excess += exc_micro[bits];
        if(pos_a->i_t + 19 <= a->t.n) {
          pos_a->leaves += count(bv_get_int(&(a->t), pos_a->i_t, 19) |
                                 (-1ULL << 19));
        } else {
          pos_a->leaves += count(bv_get_int(&(a->t), pos_a->i_t, (a->t.n - pos_a->i_t)) |
                                 (-1ULL << (a->t.n - pos_a->i_t)));
        }
      }

      // is in the last bits of the tree
      for(;;) {
        if(pos_a->i_t + 4 <= a->t.n && bv_get_int(&(a->t), pos_a->i_t, 4) == LEAF_1) {
          pos_a->leaves++;
        }
        if(bv_i(&(a->t), pos_a->i_t)) {
          pos_a->excess++;
          bv_pb(&(c->t), 1);
        } else {
          pos_a->excess--;
          bv_pb(&(c->t), 0);
        }
        pos_a->i_t++;
        if(obj_excess == pos_a->excess + 1) {
          pos_a->size = (pos_a->i_t - curr_pos) / 2;
          for(size_t i = 0; i < pos_a->leaves; i++) {
            k2bp_write_leaf(c, k2bp_read_leaf(a, pos_a->i_l + i));
            c->m += __builtin_popcount(k2bp_read_leaf(a, pos_a->i_l + i));
          }
          pos_a->i_l += pos_a->leaves;
          return;
        }
      }
    }
    for(size_t i = pos_a->i_t; i < pos_a->i_t + BLOCK_SIZE; i += 64) {
      bv_append_int(&(c->t), bv_get_int(&(a->t), pos_a->i_t, 64));
    }
    pos_a->excess = a->exc_samples[block];
    pos_a->leaves += a->leaves_samples[block];
    block++;
  }
}

static void k2bp_excdfs(k2bp_traversal_t* pos_a, const k2bp_t* a) {
  assert(bv_i(&(a->t), pos_a->i_t) == 1);
  if(bv_get_int(&(a->t), pos_a->i_t, 4) == LEAF_1) {
    pos_a->i_t += 4;
    pos_a->size = 2;
    pos_a->leaves = 1;
    return;
  }
  if(bv_get_int(&(a->t), pos_a->i_t, 2) == LEAF_0) {
    pos_a->i_t += 2;
    pos_a->size = 1;
    pos_a->leaves = 0;
    return;
  }

  pos_a->leaves = 0;
  pos_a->size = 0;
  size_t curr_pos = pos_a->i_t;
  int64_t obj_excess = pos_a->excess;
  pos_a->i_t++;
  size_t obj_pos = BLOCK_SIZE * ((pos_a->i_t + BLOCK_SIZE - 1) / BLOCK_SIZE);
  // micro tables
  for(; pos_a->i_t + 16 < obj_pos; pos_a->i_t += 16) {
    uint64_t bits = bv_get_int(&(a->t), pos_a->i_t, 16);
    // found micro block
    if(obj_excess >= pos_a->excess + exc_min_micro[bits] + 1) {
      for(uint8_t i = 0; i < 16; i++) {
        pos_a->i_t++;
        if(bits & (1ULL << i)) pos_a->excess++;
        else pos_a->excess--;

        if(obj_excess == pos_a->excess + 1) {
          pos_a->size = (pos_a->i_t - curr_pos) / 2;
          pos_a->leaves += count(bits | (-1ULL << (i + 1)));
          return;
        }
      }
    }
    pos_a->excess += exc_micro[bits];
    pos_a->leaves += count(bv_get_int(&(a->t), pos_a->i_t, 19) |
                           (-1ULL << 19));
  }
  if(pos_a->i_t % BLOCK_SIZE != 0) {
    size_t extra = BLOCK_SIZE - pos_a->i_t % BLOCK_SIZE;
    uint64_t bits = bv_get_int(&(a->t), pos_a->i_t, extra);
    size_t save_pos = pos_a->i_t;
    for(uint8_t i = 0; i < extra; i++) {
      pos_a->i_t++;
      if(bits & (1ULL << i)) pos_a->excess++;
      else pos_a->excess--;
      if(obj_excess == pos_a->excess + 1) {
        pos_a->size = (pos_a->i_t - curr_pos) / 2;
        pos_a->leaves += count(bits | (-1ULL << (i + 1)));
        return;
      }
    }
    pos_a->leaves += count(bv_get_int(&(a->t), save_pos, extra + 3) | (-1ULL << (extra + 3)));
  }

  assert(pos_a->i_t % BLOCK_SIZE == 0);

  // jumps between blocks
  size_t block = pos_a->i_t / BLOCK_SIZE;
  for(;; pos_a->i_t += BLOCK_SIZE) {
    // found block
    if(obj_excess > a->exc_min_samples[block]) {
      for(; pos_a->i_t + 16 <= a->t.n; pos_a->i_t += 16) {
        uint64_t bits = bv_get_int(&(a->t), pos_a->i_t, 16);
        // found microblock
        if(obj_excess >= pos_a->excess + exc_min_micro[bits] + 1) {
          for(uint8_t i = 0; i < 16; i++) {
            pos_a->i_t++;
            if(bits & (1ULL << i)) pos_a->excess++;
            else pos_a->excess--;
            if(obj_excess == pos_a->excess + 1) {
              pos_a->size = (pos_a->i_t - curr_pos) / 2;
              pos_a->leaves += count(bits | (-1ULL << (i + 1)));
              return;
            }
          }
        }
        pos_a->excess += exc_micro[bits];
        if(pos_a->i_t + 19 <= a->t.n) {
          pos_a->leaves += count(bv_get_int(&(a->t), pos_a->i_t, 19) |
                                 (-1ULL << 19));
        } else {
          pos_a->leaves += count(bv_get_int(&(a->t), pos_a->i_t, (a->t.n - pos_a->i_t)) |
                                 (-1ULL << (a->t.n - pos_a->i_t)));
        }
      }

      // is in the last bits of the tree
      for(;;) {
        if(pos_a->i_t + 4 <= a->t.n && bv_get_int(&(a->t), pos_a->i_t, 4) == LEAF_1) {
          pos_a->leaves++;
        }
        if(bv_i(&(a->t), pos_a->i_t)) pos_a->excess++;
        else pos_a->excess--;
        pos_a->i_t++;
        if(obj_excess == pos_a->excess + 1) {
          pos_a->size = (pos_a->i_t - curr_pos) / 2;
          return;
        }
      }
    }
    pos_a->excess = a->exc_samples[block];
    pos_a->leaves += a->leaves_samples[block];
    block++;
  }
}

static void reck2bp_mul(k2bp_traversal_t* pos_a, const k2bp_t* a, k2bp_traversal_t* pos_b, const k2bp_t* b, k2bp_t* c) {
  assert(bv_i(&(a->t), pos_a->i_t) == 1);
  assert(bv_i(&(b->t), pos_b->i_t) == 1);

  if(bv_get_int(&(a->t), pos_a->i_t, 2) == LEAF_0) { // a is full of zeros
    bv_pb(&(c->t), 1);
    bv_pb(&(c->t), 0);

    pos_a->i_t += 2;
    pos_b->i_t += pos_b->size * 2;
    pos_b->i_l += pos_b->leaves;
    return;
  }

  if(bv_get_int(&(b->t), pos_b->i_t, 2) == LEAF_0) { // a is full of zeros
    bv_pb(&(c->t), 1);
    bv_pb(&(c->t), 0);

    pos_b->i_t += 2;
    pos_a->i_t += pos_a->size * 2;
    pos_a->i_l += pos_a->leaves;
    return;
  }

  if(pos_a->msize == _K_) { // a is full of zeros
    uint8_t res = table_mul[k2bp_read_leaf(a, pos_a->i_l)][k2bp_read_leaf(b, pos_b->i_l)];
    if(res == 0) {
      bv_pb(&(c->t), 1);
      bv_pb(&(c->t), 0);
    } else {
      bv_pb(&(c->t), 1);
      bv_pb(&(c->t), 1);
      bv_pb(&(c->t), 0);
      bv_pb(&(c->t), 0);
      k2bp_write_leaf(c, res);
    }

    pos_b->i_t += 4;
    pos_a->i_t += 4;

    pos_b->i_l += 1;
    pos_a->i_l += 1;
    c->m += __builtin_popcount(res); 

    return;
  }

  if(pos_a->size < a->threshold || pos_b->size < b->threshold) {
    reck2bp_scanmul(pos_a, a, pos_b, b, c);
    return;
  }

  k2bp_traversal_t as1[4], bs1[4];
  k2bp_split(pos_a, a, as1);
  k2bp_split(pos_b, b, bs1);
  k2bp_traversal_t as2[4], bs2[4];
  k2bp_split(pos_a, a, as2);
  k2bp_split(pos_b, b, bs2);

  k2bp_t aux_c[2] = {K2BP_INITIALIZER, K2BP_INITIALIZER};
  k2bp_traversal_t aux_c_pos[2] = {K2BP_TRAVERSAL_INITIALIZER, K2BP_TRAVERSAL_INITIALIZER};
  //k2bp_t c_[4] = {K2BP_INITIALIZER, K2BP_INITIALIZER, K2BP_INITIALIZER, K2BP_INITIALIZER};
  bv_init(&(aux_c[0].t));
  aux_c[0].maxn_l = 10;
  aux_c[0].n_l = 0;
  aux_c[0].l = (uint8_t*) malloc(sizeof(uint8_t) * aux_c[0].maxn_l);
  bv_init(&(aux_c[1].t));
  aux_c[1].maxn_l = 10;
  aux_c[1].n_l = 0;
  aux_c[1].l = (uint8_t*) malloc(sizeof(uint8_t) * aux_c[1].maxn_l);

  bv_pb(&(c->t), 1);
  size_t prev_n = c->t.n;
  size_t amount_add = 0;

  reck2bp_mul(&(as1[0]), a, &(bs1[0]), b, &(aux_c[0]));
  reck2bp_mul(&(as1[1]), a, &(bs1[2]), b, &(aux_c[1]));
  aux_c[0].threshold = aux_c[0].t.n;
  aux_c[1].threshold = aux_c[1].t.n;
  //k2bp_scansum(&(aux_c[0]), &(aux_c[1]), &(c_[0]));
  k2bp_init_traversalinfo(pos_a, &(aux_c_pos[0]));
  k2bp_init_traversalinfo(pos_a, &(aux_c_pos[1]));
  reck2bp_scansum(&(aux_c_pos[0]), &(aux_c[0]), &(aux_c_pos[1]), &(aux_c[1]), c);

  aux_c[0].n_l = 0; aux_c[0].t.n = 0;
  aux_c[1].n_l = 0; aux_c[1].t.n = 0;

  reck2bp_mul(&(as2[0]), a, &(bs1[1]), b, &(aux_c[0]));
  reck2bp_mul(&(as2[1]), a, &(bs1[3]), b, &(aux_c[1]));
  aux_c[0].threshold = aux_c[0].t.n;
  aux_c[1].threshold = aux_c[1].t.n;
  //k2bp_scansum(&(aux_c[0]), &(aux_c[1]), &(c_[1]));
  k2bp_init_traversalinfo(pos_a, &(aux_c_pos[0]));
  k2bp_init_traversalinfo(pos_a, &(aux_c_pos[1]));
  reck2bp_scansum(&(aux_c_pos[0]), &(aux_c[0]), &(aux_c_pos[1]), &(aux_c[1]), c);
  aux_c[0].n_l = 0; aux_c[0].t.n = 0;
  aux_c[1].n_l = 0; aux_c[1].t.n = 0;

  reck2bp_mul(&(as1[2]), a, &(bs2[0]), b, &(aux_c[0]));
  reck2bp_mul(&(as1[3]), a, &(bs2[2]), b, &(aux_c[1]));
  aux_c[0].threshold = aux_c[0].t.n;
  aux_c[1].threshold = aux_c[1].t.n;
  //k2bp_scansum(&(aux_c[0]), &(aux_c[1]), &(c_[2]));
  k2bp_init_traversalinfo(pos_a, &(aux_c_pos[0]));
  k2bp_init_traversalinfo(pos_a, &(aux_c_pos[1]));
  reck2bp_scansum(&(aux_c_pos[0]), &(aux_c[0]), &(aux_c_pos[1]), &(aux_c[1]), c);
  aux_c[0].n_l = 0; aux_c[0].t.n = 0;
  aux_c[1].n_l = 0; aux_c[1].t.n = 0;

  reck2bp_mul(&(as2[2]), a, &(bs2[1]), b, &(aux_c[0]));
  reck2bp_mul(&(as2[3]), a, &(bs2[3]), b, &(aux_c[1]));
  aux_c[0].threshold = aux_c[0].t.n;
  aux_c[1].threshold = aux_c[1].t.n;
  //k2bp_scansum(&(aux_c[0]), &(aux_c[1]), &(c_[3]));
  k2bp_init_traversalinfo(pos_a, &(aux_c_pos[0]));
  k2bp_init_traversalinfo(pos_a, &(aux_c_pos[1]));
  reck2bp_scansum(&(aux_c_pos[0]), &(aux_c[0]), &(aux_c_pos[1]), &(aux_c[1]), c);

  k2bp_free(&(aux_c[0]));
  k2bp_free(&(aux_c[1]));
  amount_add = c->t.n - prev_n + 1;

  assert(amount_add >= 8);
  if(amount_add == 8) {
    c->t.n -= 8;
  }
  pos_a->i_t = as1[3].i_t + 1;
  pos_b->i_t = bs1[3].i_t + 1;
  pos_a->i_l = as1[3].i_l;
  pos_b->i_l = bs1[3].i_l;
  bv_pb(&(c->t), 0);
}

static void reck2bp_scanmul(k2bp_traversal_t* pos_a, const k2bp_t* a, k2bp_traversal_t* pos_b, const k2bp_t* b, k2bp_t* c) {
  assert(bv_i(&(a->t), pos_a->i_t) == 1);
  assert(bv_i(&(b->t), pos_b->i_t) == 1);

  if(bv_get_int(&(a->t), pos_a->i_t, 2) == 1) { // a is full of zeros
    bv_pb(&(c->t), 1);
    bv_pb(&(c->t), 0);

    pos_a->i_t += 2;
    if(pos_b->size > 0) {
      pos_b->i_t += pos_b->size * 2;
      pos_b->i_l += pos_b->leaves;
    } else if(b->exc_min_samples != NULL) {
      k2bp_excdfs(pos_b, b);
      pos_b->i_l += pos_b->leaves;
    } else {
      k2bp_scandfs(pos_b, b);
      pos_b->i_l += pos_b->leaves;
    }
    return;
  }

  if(bv_get_int(&(b->t), pos_b->i_t, 2) == 1) { // a is full of zeros
    bv_pb(&(c->t), 1);
    bv_pb(&(c->t), 0);

    pos_b->i_t += 2;
    if(pos_a->size > 0) {
      pos_a->i_t += pos_a->size * 2;
      pos_a->i_l += pos_a->leaves;
    } else if(a->exc_min_samples != NULL) {
      k2bp_excdfs(pos_a, a);
      pos_a->i_l += pos_a->leaves;
    } else {
      k2bp_scandfs(pos_a, a);
      pos_a->i_l += pos_a->leaves;
    }
    return;
  }

  if(pos_a->msize == _K_) { // a is full of zeros
    uint8_t res = table_mul[k2bp_read_leaf(a, pos_a->i_l)][k2bp_read_leaf(b, pos_b->i_l)];
    if(res == 0) {
      bv_pb(&(c->t), 1);
      bv_pb(&(c->t), 0);
    } else {
      bv_pb(&(c->t), 1);
      bv_pb(&(c->t), 1);
      bv_pb(&(c->t), 0);
      bv_pb(&(c->t), 0);
      k2bp_write_leaf(c, res);
      c->m += __builtin_popcount(res); 
    }

    pos_b->i_t += 4;
    pos_a->i_t += 4;

    pos_b->i_l += 1;
    pos_a->i_l += 1;

    return;
  }

  k2bp_traversal_t as1[4] = {K2BP_TRAVERSAL_INITIALIZER, K2BP_TRAVERSAL_INITIALIZER, K2BP_TRAVERSAL_INITIALIZER, K2BP_TRAVERSAL_INITIALIZER};
  k2bp_traversal_t bs1[4] = {K2BP_TRAVERSAL_INITIALIZER, K2BP_TRAVERSAL_INITIALIZER, K2BP_TRAVERSAL_INITIALIZER, K2BP_TRAVERSAL_INITIALIZER};
  k2bp_splitinfo(pos_a, a, as1);
  k2bp_splitinfo(pos_b, b, bs1);
  k2bp_traversal_t as2[4] = {K2BP_TRAVERSAL_INITIALIZER, K2BP_TRAVERSAL_INITIALIZER, K2BP_TRAVERSAL_INITIALIZER, K2BP_TRAVERSAL_INITIALIZER};
  k2bp_traversal_t bs2[4] = {K2BP_TRAVERSAL_INITIALIZER, K2BP_TRAVERSAL_INITIALIZER, K2BP_TRAVERSAL_INITIALIZER, K2BP_TRAVERSAL_INITIALIZER};
  k2bp_splitinfo(pos_a, a, as2);
  k2bp_splitinfo(pos_b, b, bs2);

  k2bp_t aux_c[3] = {K2BP_INITIALIZER, K2BP_INITIALIZER, K2BP_INITIALIZER};
  k2bp_traversal_t aux_c_pos[3] = {K2BP_TRAVERSAL_INITIALIZER, K2BP_TRAVERSAL_INITIALIZER, K2BP_TRAVERSAL_INITIALIZER};
  bv_init(&(aux_c[0].t));
  aux_c[0].maxn_l = 10;
  aux_c[0].n_l = 0;
  aux_c[0].l = (uint8_t*) malloc(sizeof(uint8_t) * aux_c[0].maxn_l);
  bv_init(&(aux_c[1].t));
  aux_c[1].maxn_l = 10;
  aux_c[1].n_l = 0;
  aux_c[1].l = (uint8_t*) malloc(sizeof(uint8_t) * aux_c[1].maxn_l);
  bv_init(&(aux_c[2].t));
  aux_c[2].maxn_l = 10;
  aux_c[2].n_l = 0;
  aux_c[2].l = (uint8_t*) malloc(sizeof(uint8_t) * aux_c[2].maxn_l);

  bv_pb(&(c->t), 1);

  size_t prev_n = c->t.n;
  if(as1[0].size == 0) {
    as1[0].i_t = pos_a->i_t + 1; as1[0].i_l = pos_a->i_l;
    as2[0].i_t = pos_a->i_t + 1; as2[0].i_l = pos_a->i_l;
  }

  if(bs1[0].size == 0) {
    bs1[0].i_t = pos_b->i_t + 1; bs1[0].i_l = pos_b->i_l;
    bs2[0].i_t = pos_b->i_t + 1; bs2[0].i_l = pos_b->i_l;
  }
  reck2bp_scanmul(&(as1[0]), a, &(bs1[0]), b, &(aux_c[0])); // X[0] = A0*B0
  as2[0].size = as1[0].size; as2[0].leaves = as1[0].leaves;
  bs2[0].size = bs1[0].size; bs2[0].leaves = bs1[0].leaves;

  as1[1].i_t = as1[0].i_t; as1[1].i_l = as1[0].i_l;
  bs1[1].i_t = bs1[0].i_t; bs1[1].i_l = bs1[0].i_l;

  as2[1].i_t = as1[0].i_t; as2[1].i_l = as1[0].i_l;
  bs2[1].i_t = bs1[0].i_t; bs2[1].i_l = bs1[0].i_l;
  reck2bp_scanmul(&(as2[0]), a, &(bs1[1]), b, &(aux_c[1])); // X = {A0*B0, A0*B1}
  bs2[1].size = bs1[1].size; bs2[1].leaves = bs1[1].leaves;

  bs1[2].i_t = bs1[1].i_t; bs1[2].i_l = bs1[1].i_l;
  bs2[2].i_t = bs1[1].i_t; bs2[2].i_l = bs1[1].i_l;
  reck2bp_scanmul(&(as1[1]), a, &(bs1[2]), b, &(aux_c[2])); // X = {A0*B0, A0*B1, A1*B2}
  as2[1].size = as1[1].size; as2[1].leaves = as1[1].leaves;

  as1[2].i_t = as1[1].i_t; as1[2].i_l = as1[1].i_l;
  as2[2].i_t = as1[1].i_t; as2[2].i_l = as1[1].i_l;

  bs1[3].i_t = bs1[2].i_t; bs1[3].i_l = bs1[2].i_l;
  bs2[3].i_t = bs1[2].i_t; bs2[3].i_l = bs1[2].i_l;

  aux_c[0].threshold = aux_c[0].t.n;
  aux_c[2].threshold = aux_c[2].t.n;
  //k2bp_scansum(&(aux_c[0]), &(aux_c[2]), &(c_[0])); // X = {_, A0*B1, _}
  k2bp_init_traversalinfo(pos_a, &(aux_c_pos[0]));
  k2bp_init_traversalinfo(pos_a, &(aux_c_pos[2]));
  reck2bp_scansum(&(aux_c_pos[0]), &(aux_c[0]), &(aux_c_pos[2]), &(aux_c[2]), c);
  aux_c[0].n_l = 0; aux_c[0].t.n = 0;
  aux_c[2].n_l = 0; aux_c[2].t.n = 0;
  k2bp_free(&(aux_c[2]));

  reck2bp_scanmul(&(as2[1]), a, &(bs1[3]), b, &(aux_c[0])); // X = {A1*B3, A0*B1, _}
  bs2[3].size = bs1[3].size; bs2[3].leaves = bs1[3].leaves;

  aux_c[0].threshold = aux_c[0].t.n;
  aux_c[1].threshold = aux_c[1].t.n;
  //k2bp_scansum(&(aux_c[0]), &(aux_c[1]), &(c_[1])); // X = {_, _, _}
  k2bp_init_traversalinfo(pos_a, &(aux_c_pos[0]));
  k2bp_init_traversalinfo(pos_a, &(aux_c_pos[1]));
  reck2bp_scansum(&(aux_c_pos[0]), &(aux_c[0]), &(aux_c_pos[1]), &(aux_c[1]), c);
  aux_c[0].n_l = 0; aux_c[0].t.n = 0;
  aux_c[1].n_l = 0; aux_c[1].t.n = 0;

  reck2bp_scanmul(&(as1[2]), a, &(bs2[0]), b, &(aux_c[0])); // X = {A2*B0, _, _}
  as2[2].size = as1[2].size; as2[2].leaves = as1[2].size;

  as1[3].i_t = as1[2].i_t; as1[3].i_l = as1[2].i_l;
  as2[3].i_t = as1[2].i_t; as2[3].i_l = as1[2].i_l;

  reck2bp_scanmul(&(as1[3]), a, &(bs2[2]), b, &(aux_c[1])); // X = {A2*B0, A3*B2, _}
  as2[3].size = as1[3].size; as2[3].leaves = as1[3].leaves;

  aux_c[0].threshold = aux_c[0].t.n;
  aux_c[1].threshold = aux_c[1].t.n;
  //k2bp_scansum(&(aux_c[0]), &(aux_c[1]), &(c_[2]));
  k2bp_init_traversalinfo(pos_a, &(aux_c_pos[0]));
  k2bp_init_traversalinfo(pos_a, &(aux_c_pos[1]));
  reck2bp_scansum(&(aux_c_pos[0]), &(aux_c[0]), &(aux_c_pos[1]), &(aux_c[1]), c);
  aux_c[0].n_l = 0; aux_c[0].t.n = 0;
  aux_c[1].n_l = 0; aux_c[1].t.n = 0;

  reck2bp_scanmul(&(as2[2]), a, &(bs2[1]), b, &(aux_c[0]));
  reck2bp_scanmul(&(as2[3]), a, &(bs2[3]), b, &(aux_c[1]));
  aux_c[0].threshold = aux_c[0].t.n;
  aux_c[1].threshold = aux_c[1].t.n;
  //k2bp_scansum(&(aux_c[0]), &(aux_c[1]), &(c_[3]));
  k2bp_init_traversalinfo(pos_a, &(aux_c_pos[0]));
  k2bp_init_traversalinfo(pos_a, &(aux_c_pos[1]));
  reck2bp_scansum(&(aux_c_pos[0]), &(aux_c[0]), &(aux_c_pos[1]), &(aux_c[1]), c);

  k2bp_free(&(aux_c[0]));
  k2bp_free(&(aux_c[1]));

  size_t amount_add = c->t.n - prev_n + 1;
  assert(amount_add >= 8);
  if(amount_add == 8) {
    c->t.n -= 8;
  }
  pos_a->i_t = as1[3].i_t + 1;
  pos_b->i_t = bs1[3].i_t + 1;
  pos_a->i_l = as1[3].i_l;
  pos_b->i_l = bs1[3].i_l;
  bv_pb(&(c->t), 0);
}

static void reck2bp_scansum(k2bp_traversal_t* pos_a, const k2bp_t* a, k2bp_traversal_t* pos_b, const k2bp_t* b, k2bp_t* c) {
  assert(bv_i(&(a->t), pos_a->i_t) == 1);
  assert(bv_i(&(b->t), pos_b->i_t) == 1);
  if(bv_get_int(&(a->t), pos_a->i_t, 2) == LEAF_0) {
    if(pos_b->size > 0) {
      size_t end_tree = pos_b->i_t + pos_b->size * 2;
      for(; pos_b->i_t + 64 < end_tree; pos_b->i_t += 64) {
        bv_append_int(&(c->t), bv_get_int(&(b->t), pos_b->i_t, 64));
        assert(bv_get_int(&(c->t), c->t.n - 64, 64) == bv_get_int(&(b->t), pos_b->i_t, 64));
      }
      for(;pos_b->i_t < end_tree; pos_b->i_t++) {
        bv_pb(&(c->t), bv_i(&(b->t), pos_b->i_t));
        assert(bv_i(&(c->t), c->t.n - 1) == bv_i(&(b->t), pos_b->i_t));
      }
      for(size_t i = 0; i < pos_b->leaves; i++) {
        k2bp_write_leaf(c, k2bp_read_leaf(b, pos_b->i_l + i));
        assert(k2bp_read_leaf(c, c->n_l - 1) == k2bp_read_leaf(b, pos_b->i_l + i));
        c->m += __builtin_popcount(k2bp_read_leaf(c, c->n_l - 1));
      }
      pos_b->i_l += pos_b->leaves;
    } else {
      k2bp_scandfs_copy(pos_b, b, c);
    }
    pos_a->i_t += 2;
    return;
  }

  if(bv_get_int(&(b->t), pos_b->i_t, 2) == LEAF_0) {
    if(pos_a->size > 0) {
      size_t end_tree = pos_a->i_t + pos_a->size * 2;
      for(; pos_a->i_t + 64 < end_tree; pos_a->i_t += 64) {
        bv_append_int(&(c->t), bv_get_int(&(a->t), pos_a->i_t, 64));
        assert(bv_get_int(&(c->t), c->t.n - 64, 64) == bv_get_int(&(a->t), pos_a->i_t, 64));
      }
      for(;pos_a->i_t < end_tree; pos_a->i_t++) {
        bv_pb(&(c->t), bv_i(&(a->t), pos_a->i_t));
        assert(bv_i(&(c->t), c->t.n - 1) == bv_i(&(a->t), pos_a->i_t));
      }
      for(size_t i = 0; i < pos_a->leaves; i++) {
        k2bp_write_leaf(c, k2bp_read_leaf(a, pos_a->i_l + i));
        assert(k2bp_read_leaf(c, c->n_l - 1) == k2bp_read_leaf(a, pos_a->i_l + i));
        c->m += __builtin_popcount(k2bp_read_leaf(c, c->n_l - 1));
      }
      pos_a->i_l += pos_a->leaves;
    } else {
      k2bp_scandfs_copy(pos_a, a, c);
    }
    pos_b->i_t += 2;
    return;
  }


  if(bv_get_int(&(a->t), pos_a->i_t, 4) == LEAF_1 && bv_get_int(&(b->t), pos_b->i_t, 4) == LEAF_1) {
    bv_pb(&(c->t), 1);
    bv_pb(&(c->t), 1);
    bv_pb(&(c->t), 0);
    bv_pb(&(c->t), 0);

    k2bp_write_leaf(c, k2bp_read_leaf(a, pos_a->i_l) | k2bp_read_leaf(b, pos_b->i_l));
    pos_a->i_t += 4;
    pos_b->i_t += 4;
    pos_b->i_l += 1;
    pos_a->i_l += 1;
    c->m += __builtin_popcount(k2bp_read_leaf(c, c->n_l - 1));
    return;
  }

  k2bp_traversal_t aux_a[4] = {K2BP_TRAVERSAL_INITIALIZER,
                               K2BP_TRAVERSAL_INITIALIZER,
                               K2BP_TRAVERSAL_INITIALIZER,
                               K2BP_TRAVERSAL_INITIALIZER};
  k2bp_traversal_t aux_b[4] = {K2BP_TRAVERSAL_INITIALIZER,
                               K2BP_TRAVERSAL_INITIALIZER,
                               K2BP_TRAVERSAL_INITIALIZER,
                               K2BP_TRAVERSAL_INITIALIZER};

  k2bp_splitinfo(pos_a, a, aux_a);
  k2bp_splitinfo(pos_b, b, aux_b);

  aux_a[0].i_t = pos_a->i_t + 1;
  aux_a[0].i_l = pos_a->i_l;

  aux_b[0].i_t = pos_b->i_t + 1;
  aux_b[0].i_l = pos_b->i_l;

  bv_pb(&(c->t), 1);
  reck2bp_scansum(&(aux_a[0]), a, &(aux_b[0]), b, c);
  aux_a[1].i_t = aux_a[0].i_t;
  aux_a[1].i_l = aux_a[0].i_l;

  aux_b[1].i_t = aux_b[0].i_t;
  aux_b[1].i_l = aux_b[0].i_l;

  reck2bp_scansum(&(aux_a[1]), a, &(aux_b[1]), b, c);
  aux_a[2].i_t = aux_a[1].i_t;
  aux_a[2].i_l = aux_a[1].i_l;

  aux_b[2].i_t = aux_b[1].i_t;
  aux_b[2].i_l = aux_b[1].i_l;

  reck2bp_scansum(&(aux_a[2]), a, &(aux_b[2]), b, c);
  aux_a[3].i_t = aux_a[2].i_t;
  aux_a[3].i_l = aux_a[2].i_l;

  aux_b[3].i_t = aux_b[2].i_t;
  aux_b[3].i_l = aux_b[2].i_l;
  reck2bp_scansum(&(aux_a[3]), a, &(aux_b[3]), b, c);

  pos_a->i_t = aux_a[3].i_t + 1;
  pos_a->i_l = aux_a[3].i_l;
  
  pos_b->i_t = aux_b[3].i_t + 1;
  pos_b->i_l = aux_b[3].i_l;
  bv_pb(&(c->t), 0);
}


static void reck2bp_sum(k2bp_traversal_t* pos_a, const k2bp_t* a, k2bp_traversal_t* pos_b, const k2bp_t* b, k2bp_t* c) {

  assert(bv_i(&(a->t), pos_a->i_t) == 1);
  assert(bv_i(&(b->t), pos_b->i_t) == 1);
  assert(pos_a->size > 0);
  assert(pos_b->size > 0);
  if(bv_get_int(&(a->t), pos_a->i_t, 2) == LEAF_0) {
    size_t end_tree = pos_b->i_t + pos_b->size * 2;
    for(; pos_b->i_t + 64 < end_tree; pos_b->i_t += 64) {
      bv_append_int(&(c->t), bv_get_int(&(b->t), pos_b->i_t, 64));
      assert(bv_get_int(&(c->t), c->t.n - 64, 64) == bv_get_int(&(b->t), pos_b->i_t, 64));
    }
    for(;pos_b->i_t < end_tree; pos_b->i_t++) {
      bv_pb(&(c->t), bv_i(&(b->t), pos_b->i_t));
      assert(bv_i(&(c->t), c->t.n - 1) == bv_i(&(b->t), pos_b->i_t));
    }
    for(size_t i = 0; i < pos_b->leaves; i++) {
      k2bp_write_leaf(c, k2bp_read_leaf(b, pos_b->i_l + i));
      assert(k2bp_read_leaf(c, c->n_l - 1) == k2bp_read_leaf(b, pos_b->i_l + i));
      c->m += __builtin_popcount(k2bp_read_leaf(c, c->n_l - 1));
    }
    pos_a->i_t += 2;
    pos_b->i_l += pos_b->leaves;
    return;
  }

  if(bv_get_int(&(b->t), pos_b->i_t, 2) == LEAF_0) {
    size_t end_tree = pos_a->i_t + pos_a->size * 2;
    for(; pos_a->i_t + 64 < end_tree; pos_a->i_t += 64) {
      bv_append_int(&(c->t), bv_get_int(&(a->t), pos_a->i_t, 64));
      assert(bv_get_int(&(c->t), c->t.n - 64, 64) == bv_get_int(&(a->t), pos_a->i_t, 64));
    }
    for(;pos_a->i_t < end_tree; pos_a->i_t++) {
      bv_pb(&(c->t), bv_i(&(a->t), pos_a->i_t));
      assert(bv_i(&(c->t), c->t.n - 1) == bv_i(&(a->t), pos_a->i_t));
    }
    for(size_t i = 0; i < pos_a->leaves; i++) {
      k2bp_write_leaf(c, k2bp_read_leaf(a, pos_a->i_l + i));
      assert(k2bp_read_leaf(c, c->n_l - 1) == k2bp_read_leaf(a, pos_a->i_l + i));
      c->m += __builtin_popcount(k2bp_read_leaf(c, c->n_l - 1));
    }
    pos_b->i_t += 2;
    pos_a->i_l += pos_a->leaves;
    return;
  }


  if(bv_get_int(&(a->t), pos_a->i_t, 4) == LEAF_1 && bv_get_int(&(b->t), pos_b->i_t, 4) == LEAF_1) {
    bv_pb(&(c->t), 1);
    bv_pb(&(c->t), 1);
    bv_pb(&(c->t), 0);
    bv_pb(&(c->t), 0);

    k2bp_write_leaf(c, k2bp_read_leaf(a, pos_a->i_l) | k2bp_read_leaf(b, pos_b->i_l));
    pos_a->i_t += 4;
    pos_b->i_t += 4;
    pos_b->i_l += 1;
    pos_a->i_l += 1;
    c->m += __builtin_popcount(k2bp_read_leaf(c, c->n_l - 1));
    return;
  }

  k2bp_traversal_t aux_a[4] = {K2BP_TRAVERSAL_INITIALIZER,
                               K2BP_TRAVERSAL_INITIALIZER,
                               K2BP_TRAVERSAL_INITIALIZER,
                               K2BP_TRAVERSAL_INITIALIZER};
  k2bp_traversal_t aux_b[4] = {K2BP_TRAVERSAL_INITIALIZER,
                               K2BP_TRAVERSAL_INITIALIZER,
                               K2BP_TRAVERSAL_INITIALIZER,
                               K2BP_TRAVERSAL_INITIALIZER};
  k2bp_split(pos_a, a, aux_a);
  k2bp_split(pos_b, b, aux_b);

  bv_pb(&(c->t), 1);
  reck2bp_sum(&(aux_a[0]), a, &(aux_b[0]), b, c);
  assert(aux_a[0].i_l == aux_a[1].i_l);
  assert(aux_b[0].i_l == aux_b[1].i_l);
  reck2bp_sum(&(aux_a[1]), a, &(aux_b[1]), b, c);
  assert(aux_a[1].i_l == aux_a[2].i_l);
  assert(aux_b[1].i_l == aux_b[2].i_l);
  reck2bp_sum(&(aux_a[2]), a, &(aux_b[2]), b, c);
  assert(aux_a[2].i_l == aux_a[3].i_l);
  assert(aux_b[2].i_l == aux_b[3].i_l);
  reck2bp_sum(&(aux_a[3]), a, &(aux_b[3]), b, c);
  pos_a->i_t = aux_a[3].i_t + 1;
  pos_a->i_l = aux_a[3].i_l;
  
  pos_b->i_t = aux_b[3].i_t + 1;
  pos_b->i_l = aux_b[3].i_l;
  bv_pb(&(c->t), 0);
}

static void k2bp_init_traversalinfo(const k2bp_traversal_t* pos_a, k2bp_traversal_t* copy_a) {
  if(pos_a->msize != 0) { // used
    copy_a->msize = pos_a->msize;
    copy_a->x = pos_a->x;
    copy_a->y = pos_a->y;
  }
  copy_a->excess = 1;
  copy_a->i_t = 0;
  copy_a->i_l = 0;
  copy_a->size = 0;
  copy_a->node = 0;
  copy_a->leaves = 0;
  copy_a->i_p = 0;
  copy_a->rank_pointers = pos_a->rank_pointers;
}

static void k2bp_splitinfo(const k2bp_traversal_t* pos_a, const k2bp_t* a, k2bp_traversal_t* splits) {
  if(pos_a->msize != 0) { // used
    splits[0].msize = pos_a->msize / 2;
    splits[0].x = pos_a->x;
    splits[0].y = pos_a->y;

    splits[1].msize = pos_a->msize / 2;
    splits[1].x = pos_a->x;
    splits[1].y = pos_a->y + pos_a->msize / 2;

    splits[2].msize = pos_a->msize / 2;
    splits[2].x = pos_a->x + pos_a->msize / 2;
    splits[2].y = pos_a->y;

    splits[3].msize = pos_a->msize / 2;
    splits[3].x = pos_a->x + pos_a->msize / 2;
    splits[3].y = pos_a->y + pos_a->msize / 2;
  }

  splits[0].excess = pos_a->excess + 1;
  splits[1].excess = pos_a->excess + 1;
  splits[2].excess = pos_a->excess + 1;
  splits[3].excess = pos_a->excess + 1;

  if(a->subtreeinfo != NULL && pos_a->size >= a->threshold) {
    splits[0].i_t = pos_a->i_t + 1;
    splits[0].i_l = pos_a->i_l;
    splits[0].node = pos_a->node + 3;
    splits[0].size = GET_NODES(a->subtreeinfo[pos_a->node]);
    splits[0].leaves = a->leavesinfo[pos_a->node];

    splits[1].i_t = splits[0].i_t + splits[0].size * 2;
    splits[1].i_l = splits[0].i_l + splits[0].leaves;
    splits[1].node = splits[0].node + GET_SKIPS(a->subtreeinfo[pos_a->node]);
    splits[1].size = GET_NODES(a->subtreeinfo[pos_a->node + 1]);
    splits[1].leaves = a->leavesinfo[pos_a->node + 1];

    splits[2].i_t = splits[1].i_t + splits[1].size * 2;
    splits[2].i_l = splits[1].i_l + splits[1].leaves;
    splits[2].node = splits[1].node + GET_SKIPS(a->subtreeinfo[pos_a->node + 1]);
    splits[2].size = GET_NODES(a->subtreeinfo[pos_a->node + 2]);
    splits[2].leaves = a->leavesinfo[pos_a->node + 2];

    splits[3].i_t = splits[2].i_t + splits[2].size * 2;
    splits[3].i_l = splits[2].i_l + splits[2].leaves;
    splits[3].node = splits[2].node + GET_SKIPS(a->subtreeinfo[pos_a->node + 2]);
    splits[3].size = pos_a->size -
                     (splits[0].size + splits[1].size + splits[2].size + 1);
    splits[3].leaves = pos_a->leaves -
                       (splits[0].leaves + splits[1].leaves + splits[2].leaves);
  }
}

static void k2bp_split(const k2bp_traversal_t* pos_a, const k2bp_t* a, k2bp_traversal_t* splits) {
  if(pos_a->msize != 0) { // used
    splits[0].msize = pos_a->msize / 2;
    splits[0].x = pos_a->x;
    splits[0].y = pos_a->y;

    splits[1].msize = pos_a->msize / 2;
    splits[1].x = pos_a->x;
    splits[1].y = pos_a->y + pos_a->msize / 2;

    splits[2].msize = pos_a->msize / 2;
    splits[2].x = pos_a->x + pos_a->msize / 2;
    splits[2].y = pos_a->y;

    splits[3].msize = pos_a->msize / 2;
    splits[3].x = pos_a->x + pos_a->msize / 2;
    splits[3].y = pos_a->y + pos_a->msize / 2;
  }

  splits[0].excess = pos_a->excess + 1;
  splits[1].excess = pos_a->excess + 1;
  splits[2].excess = pos_a->excess + 1;
  splits[3].excess = pos_a->excess + 1;

  if(a->subtreeinfo != NULL && pos_a->size >= a->threshold) {
    splits[0].i_t = pos_a->i_t + 1;
    splits[0].i_l = pos_a->i_l;
    splits[0].node = pos_a->node + 3;
    splits[0].size = GET_NODES(a->subtreeinfo[pos_a->node]);
    splits[0].leaves = a->leavesinfo[pos_a->node];

    splits[1].i_t = splits[0].i_t + splits[0].size * 2;
    splits[1].i_l = splits[0].i_l + splits[0].leaves;
    splits[1].node = splits[0].node + GET_SKIPS(a->subtreeinfo[pos_a->node]);
    splits[1].size = GET_NODES(a->subtreeinfo[pos_a->node + 1]);
    splits[1].leaves = a->leavesinfo[pos_a->node + 1];

    splits[2].i_t = splits[1].i_t + splits[1].size * 2;
    splits[2].i_l = splits[1].i_l + splits[1].leaves;
    splits[2].node = splits[1].node + GET_SKIPS(a->subtreeinfo[pos_a->node + 1]);
    splits[2].size = GET_NODES(a->subtreeinfo[pos_a->node + 2]);
    splits[2].leaves = a->leavesinfo[pos_a->node + 2];

    splits[3].i_t = splits[2].i_t + splits[2].size * 2;
    splits[3].i_l = splits[2].i_l + splits[2].leaves;
    splits[3].node = splits[2].node + GET_SKIPS(a->subtreeinfo[pos_a->node + 2]);
    splits[3].size = pos_a->size -
                     (splits[0].size + splits[1].size + splits[2].size + 1);
    splits[3].leaves = pos_a->leaves -
                       (splits[0].leaves + splits[1].leaves + splits[2].leaves);
  } else if(pos_a->size >= BLOCK_SIZE && a->exc_min_samples != NULL) {
    k2bp_traversal_t aux_pos = K2BP_TRAVERSAL_INITIALIZER;
    aux_pos.i_t = pos_a->i_t + 1;
    aux_pos.i_l = pos_a->i_l;
    aux_pos.excess = pos_a->excess + 1;

    splits[0].i_t = aux_pos.i_t;
    splits[0].i_l = aux_pos.i_l;
    k2bp_excdfs(&aux_pos, a);
    splits[0].size = aux_pos.size;
    splits[0].leaves = aux_pos.leaves;
    aux_pos.i_l += aux_pos.leaves;
    aux_pos.excess = pos_a->excess + 1;

    splits[1].i_t = aux_pos.i_t;
    splits[1].i_l = aux_pos.i_l;
    k2bp_excdfs(&aux_pos, a);
    splits[1].size = aux_pos.size;
    splits[1].leaves = aux_pos.leaves;
    aux_pos.i_l += aux_pos.leaves;
    aux_pos.excess = pos_a->excess + 1;

    splits[2].i_t = aux_pos.i_t;
    splits[2].i_l = aux_pos.i_l;
    k2bp_excdfs(&aux_pos, a);
    splits[2].size = aux_pos.size;
    splits[2].leaves = aux_pos.leaves;
    aux_pos.i_l += aux_pos.leaves;
    aux_pos.excess = pos_a->excess + 1;

    splits[3].i_t = aux_pos.i_t;
    splits[3].i_l = aux_pos.i_l;
    k2bp_excdfs(&aux_pos, a);
    splits[3].size = aux_pos.size;
    splits[3].leaves = aux_pos.leaves;
  } else {
    k2bp_traversal_t aux_pos = K2BP_TRAVERSAL_INITIALIZER;
    aux_pos.i_t = pos_a->i_t + 1;
    aux_pos.i_l = pos_a->i_l;
    aux_pos.excess = pos_a->excess + 1;

    splits[0].i_t = aux_pos.i_t;
    splits[0].i_l = aux_pos.i_l;
    k2bp_scandfs(&aux_pos, a);
    splits[0].size = aux_pos.size;
    splits[0].leaves = aux_pos.leaves;
    aux_pos.i_l += aux_pos.leaves;
    aux_pos.excess = pos_a->excess + 1;

    splits[1].i_t = aux_pos.i_t;
    splits[1].i_l = aux_pos.i_l;
    k2bp_scandfs(&aux_pos, a);
    splits[1].size = aux_pos.size;
    splits[1].leaves = aux_pos.leaves;
    aux_pos.i_l += aux_pos.leaves;
    aux_pos.excess = pos_a->excess + 1;

    splits[2].i_t = aux_pos.i_t;
    splits[2].i_l = aux_pos.i_l;
    k2bp_scandfs(&aux_pos, a);
    splits[2].size = aux_pos.size;
    splits[2].leaves = aux_pos.leaves;
    aux_pos.i_l += aux_pos.leaves;
    aux_pos.excess = pos_a->excess + 1;

    splits[3].i_t = aux_pos.i_t;
    splits[3].i_l = aux_pos.i_l;
    k2bp_scandfs(&aux_pos, a);
    splits[3].size = aux_pos.size;
    splits[3].leaves = aux_pos.leaves;
  }
}

static void reck2bp_checksubtree_info(k2bp_traversal_t* pos_a, const k2bp_t* a) {
  assert(bv_i(&(a->t), pos_a->i_t) == 1);

  if(pos_a->i_t + 3 < a->t.n && bv_get_int(&(a->t), pos_a->i_t, 4) == LEAF_1) {
    pos_a->i_t += 4;
    pos_a->i_l++;
    return;
  }

  if(pos_a->i_t + 1 < a->t.n && bv_get_int(&(a->t), pos_a->i_t, 2) == LEAF_0) {
    pos_a->i_t += 2;
    return;
  }
  
  if(pos_a->size < a->threshold) {
    pos_a->i_t++;
    reck2bp_checksubtree_info(pos_a, a);
    reck2bp_checksubtree_info(pos_a, a);
    reck2bp_checksubtree_info(pos_a, a);
    reck2bp_checksubtree_info(pos_a, a);
    pos_a->i_t++;
    return;
  }

  uint64_t c_sizes[4] = {0, 0, 0, 0};
  uint32_t c_leaves[4] = {0, 0, 0, 0};

  pos_a->i_t++;
  k2bp_traversal_t aux_pos = {0, 0, 0, pos_a->i_t, pos_a->i_l, pos_a->node, pos_a->size, pos_a->leaves, 0};
  pos_a->node += 3;
  pos_a->size = GET_NODES(a->subtreeinfo[aux_pos.node]);
  pos_a->leaves = a->leavesinfo[aux_pos.node];
  reck2bp_checksubtree_info(pos_a, a);
  c_sizes[0] = (pos_a->i_t - aux_pos.i_t) / 2;
  c_leaves[0] = pos_a->i_l - aux_pos.i_l;

  aux_pos.i_t = pos_a->i_t;
  aux_pos.i_l = pos_a->i_l;
  pos_a->node = aux_pos.node + 3 + GET_SKIPS(a->subtreeinfo[aux_pos.node]);
  pos_a->size = GET_NODES(a->subtreeinfo[aux_pos.node + 1]);
  pos_a->leaves = a->leavesinfo[aux_pos.node + 1];
  reck2bp_checksubtree_info(pos_a, a);
  c_sizes[1] = (pos_a->i_t - aux_pos.i_t) / 2;
  c_leaves[1] = pos_a->i_l - aux_pos.i_l;

  aux_pos.i_t = pos_a->i_t;
  aux_pos.i_l = pos_a->i_l;
  pos_a->node = aux_pos.node + 3 + GET_SKIPS(a->subtreeinfo[aux_pos.node])
                                 + GET_SKIPS(a->subtreeinfo[aux_pos.node + 1]);
  pos_a->size = GET_NODES(a->subtreeinfo[aux_pos.node + 2]);
  pos_a->leaves = a->leavesinfo[aux_pos.node + 2];
  reck2bp_checksubtree_info(pos_a, a);
  c_sizes[2] = (pos_a->i_t - aux_pos.i_t) / 2;
  c_leaves[2] = pos_a->i_l - aux_pos.i_l;

  aux_pos.i_t = pos_a->i_t;
  aux_pos.i_l = pos_a->i_l;
  pos_a->node = aux_pos.node + 3 + GET_SKIPS(a->subtreeinfo[aux_pos.node])
                                 + GET_SKIPS(a->subtreeinfo[aux_pos.node + 1])
         + GET_SKIPS(a->subtreeinfo[aux_pos.node + 2]);
  pos_a->size = aux_pos.size - (GET_NODES(a->subtreeinfo[aux_pos.node])
                + GET_NODES(a->subtreeinfo[aux_pos.node + 1])
                + GET_NODES(a->subtreeinfo[aux_pos.node + 2]) + 1);
  pos_a->leaves = aux_pos.leaves - (a->leavesinfo[aux_pos.node]
         + a->leavesinfo[aux_pos.node + 1]
                  + a->leavesinfo[aux_pos.node + 2]);
  reck2bp_checksubtree_info(pos_a, a);
  c_sizes[3] = (pos_a->i_t - aux_pos.i_t) / 2;
  c_leaves[3] = pos_a->i_l - aux_pos.i_l;

  pos_a->i_t++;
  if(c_sizes[0] + c_sizes[1] + c_sizes[2] + c_sizes[3] + 1 >= a->threshold) {
    uint64_t accum_size = GET_NODES(a->subtreeinfo[aux_pos.node]);
    uint32_t leaves = a->leavesinfo[aux_pos.node];
    
    assert(c_sizes[0] == GET_NODES(a->subtreeinfo[aux_pos.node]));
    assert(c_leaves[0] == a->leavesinfo[aux_pos.node]);

    accum_size += GET_NODES(a->subtreeinfo[aux_pos.node + 1]);
    leaves += a->leavesinfo[aux_pos.node + 1];

    assert(c_sizes[1] == GET_NODES(a->subtreeinfo[aux_pos.node + 1]));
    assert(c_leaves[1] == a->leavesinfo[aux_pos.node + 1]);

    accum_size += GET_NODES(a->subtreeinfo[aux_pos.node + 2]);
    leaves += a->leavesinfo[aux_pos.node + 2];

    assert(c_sizes[2] == GET_NODES(a->subtreeinfo[aux_pos.node + 2]));
    assert(c_leaves[2] == a->leavesinfo[aux_pos.node + 2]);

    assert(c_sizes[3] == aux_pos.size - accum_size - 1);
    assert(c_leaves[3] == aux_pos.leaves - leaves);
  }
}

static void reck2bp_addsubtree_info(k2bp_traversal_t* pos_a, k2bp_t* a, uint32_t* leaves, vu64_t* subinfo, vu64_t* leavesinfo) {
  assert(bv_i(&(a->t), pos_a->i_t) == 1);

  if(pos_a->i_t + 3 < a->t.n && bv_get_int(&(a->t), pos_a->i_t, 4) == LEAF_1) {
    pos_a->i_t += 4;
    (*leaves)++;
    return;
  }

  if(pos_a->i_t + 1 < a->t.n && bv_get_int(&(a->t), pos_a->i_t, 2) == LEAF_0) {
    pos_a->i_t += 2;
    return;
  }

  vu64_t c_subtreeinfo[4];
  vu64_t c_leavesinfo[4];
  for(size_t i = 0; i < 4; i++) {
    vu64_init(&(c_subtreeinfo[i]));
    vu64_init(&(c_leavesinfo[i]));
  }

  uint64_t c_sizes[4] = {0, 0, 0, 0};
  uint32_t c_leaves[4] = {0, 0, 0, 0};

  pos_a->i_t++;
  k2bp_traversal_t aux_pos = {0, 0, 0, pos_a->i_t, 0, 0, 0, 0, 0};
  reck2bp_addsubtree_info(pos_a, a, &(c_leaves[0]), &(c_subtreeinfo[0]), &(c_leavesinfo[0]));
  c_sizes[0] = (pos_a->i_t - aux_pos.i_t) / 2;
  assert(c_subtreeinfo[0].n == c_leavesinfo[0].n);

  aux_pos.i_t = pos_a->i_t;
  reck2bp_addsubtree_info(pos_a, a, &(c_leaves[1]), &(c_subtreeinfo[1]), &(c_leavesinfo[1]));
  c_sizes[1] = (pos_a->i_t - aux_pos.i_t) / 2;
  assert(c_subtreeinfo[1].n == c_leavesinfo[1].n);

  aux_pos.i_t = pos_a->i_t;
  reck2bp_addsubtree_info(pos_a, a, &(c_leaves[2]), &(c_subtreeinfo[2]), &(c_leavesinfo[2]));
  c_sizes[2] = (pos_a->i_t - aux_pos.i_t) / 2;
  assert(c_subtreeinfo[2].n == c_leavesinfo[2].n);

  aux_pos.i_t = pos_a->i_t;
  reck2bp_addsubtree_info(pos_a, a, &(c_leaves[3]), &(c_subtreeinfo[3]), &(c_leavesinfo[3]));
  c_sizes[3] = (pos_a->i_t - aux_pos.i_t) / 2;
  assert(c_subtreeinfo[3].n == c_leavesinfo[3].n);

  pos_a->i_t++;
  if(c_sizes[0] + c_sizes[1] + c_sizes[2] + c_sizes[3] + 1 >= a->threshold) {
    size_t curr_i = subinfo->n;
    vu64_grow(subinfo, 3 + c_subtreeinfo[0].n + c_subtreeinfo[1].n +
                           c_subtreeinfo[2].n + c_subtreeinfo[3].n);
    vu64_grow(leavesinfo, 3 + c_leavesinfo[0].n + c_leavesinfo[1].n +
                              c_leavesinfo[2].n + c_leavesinfo[3].n);

    subinfo->v[curr_i] = ENCODE(c_subtreeinfo[0].n, c_sizes[0]);
    leavesinfo->v[curr_i] = c_leaves[0];
    curr_i++;

    subinfo->v[curr_i] = ENCODE(c_subtreeinfo[1].n, c_sizes[1]);
    leavesinfo->v[curr_i] = c_leaves[1];
    curr_i++;

    subinfo->v[curr_i] = ENCODE(c_subtreeinfo[2].n, c_sizes[2]);
    leavesinfo->v[curr_i] = c_leaves[2];
    curr_i++;

    if(c_sizes[0] >= a->threshold) {
      for(size_t i = 0; i < c_subtreeinfo[0].n; i++) {
        subinfo->v[curr_i] = c_subtreeinfo[0].v[i];
        leavesinfo->v[curr_i] = c_leavesinfo[0].v[i];
        curr_i++;
      }
    }

    if(c_sizes[1] >= a->threshold) {
      for(size_t i = 0; i < c_subtreeinfo[1].n; i++) {
        subinfo->v[curr_i] = c_subtreeinfo[1].v[i];
        leavesinfo->v[curr_i] = c_leavesinfo[1].v[i];
        curr_i++;
      }
    }

    if(c_sizes[2] >= a->threshold) {
      for(size_t i = 0; i < c_subtreeinfo[2].n; i++) {
        subinfo->v[curr_i] = c_subtreeinfo[2].v[i];
        leavesinfo->v[curr_i] = c_leavesinfo[2].v[i];
        curr_i++;
      }
    }

    if(c_sizes[3] >= a->threshold) {
      for(size_t i = 0; i < c_subtreeinfo[3].n; i++) {
        subinfo->v[curr_i] = c_subtreeinfo[3].v[i];
        leavesinfo->v[curr_i] = c_leavesinfo[3].v[i];
        curr_i++;
      }
    }
    assert(curr_i == subinfo->n);
  }
  *leaves = c_leaves[0] + c_leaves[1] + c_leaves[2] + c_leaves[3];
  vu64_free(&(c_subtreeinfo[0]));
  vu64_free(&(c_subtreeinfo[1]));
  vu64_free(&(c_subtreeinfo[2]));
  vu64_free(&(c_subtreeinfo[3]));

  vu64_free(&(c_leavesinfo[0]));
  vu64_free(&(c_leavesinfo[1]));
  vu64_free(&(c_leavesinfo[2]));
  vu64_free(&(c_leavesinfo[3]));
}

static void reck2bp_nonzeros(k2bp_traversal_t* pos_a, const k2bp_t* a, uint32_t* arr, size_t* n) {
  assert(pos_a->i_t < a->t.n && pos_a->i_l <= a->n_l);
  assert(bv_i(&(a->t), pos_a->i_t) == 1);
  if(pos_a->i_t + 1 < a->t.n && bv_i(&(a->t), pos_a->i_t + 1) == 0) {
    pos_a->i_t += 2;
    return;
  }
  
  if(pos_a->msize == _K_) {
    uint8_t leaf = k2bp_read_leaf(a, pos_a->i_l);
    pos_a->i_l++;
    pos_a->i_t += 4;
    if(leaf & 1) {
      arr[*n] = pos_a->x;
      arr[*n + 1] = pos_a->y;
      *n += 2;
    }
    if(leaf & 2) {
      arr[*n] = pos_a->x;
      arr[*n + 1] = pos_a->y + 1; 
      *n += 2;
    }
    if(leaf & 4) {
      arr[*n] = pos_a->x + 1;
      arr[*n + 1] = pos_a->y; 
      *n += 2;
    }
    if(leaf & 8) {
      arr[*n] = pos_a->x + 1;
      arr[*n + 1] = pos_a->y + 1; 
      *n += 2;
    }
    return;
  }

  k2bp_traversal_t pos_aux = {pos_a->msize / 2, pos_a->x, pos_a->y, pos_a->i_t + 1, pos_a->i_l, 0, 0, 0, 0};
  reck2bp_nonzeros(&pos_aux, a, arr, n);

  pos_aux.x = pos_a->x;
  pos_aux.y = pos_a->y + pos_a->msize / 2;
  reck2bp_nonzeros(&pos_aux, a, arr, n);

  pos_aux.x = pos_a->x + pos_a->msize / 2;
  pos_aux.y = pos_a->y;
  reck2bp_nonzeros(&pos_aux, a, arr, n);

  pos_aux.x = pos_a->x + pos_a->msize / 2;
  pos_aux.y = pos_a->y + pos_a->msize / 2;
  reck2bp_nonzeros(&pos_aux, a, arr, n);
  pos_a->i_t = pos_aux.i_t + 1;
  pos_a->i_l = pos_aux.i_l;
}

static uint8_t encode_leaf(uint64_t ia[], size_t n, size_t smin) {
  uint8_t leaf = 0;
  for(size_t i = 0; i < n; i++) {
    size_t pos = (ia[i] - smin);
    leaf |= 1 << pos;
  }
  return leaf;
}

static void reccreate_k2bp(uint64_t ia[], size_t n, size_t smin, size_t size, k2bp_t *a) {
  uint64_t range = (size / 2) * (size / 2);
  uint64_t left = smin + range;
  uint64_t mid = left + range;
  uint64_t right = mid + range;

  size_t imid = binsearch(ia, n, mid);
  size_t ileft = imid >0 ? binsearch(ia, imid, left) : 0;
  size_t iright = imid < n ? binsearch(ia + imid, n - imid, right) + imid : n;

  bv_pb(&(a->t), 1); // you know that has at least 1 one
  if(ileft > 0) { // a[00]
    if(size == 2 * _K_) { // leaf
      bv_pb(&(a->t), 1);
      bv_pb(&(a->t), 1);
      bv_pb(&(a->t), 0);
      bv_pb(&(a->t), 0);
      uint8_t leaf = encode_leaf(ia, ileft, smin);
      k2bp_write_leaf(a, leaf);
    } else {
      reccreate_k2bp(ia, ileft, smin, size / 2, a);
    }
  } else {
    bv_pb(&(a->t), 1);
    bv_pb(&(a->t), 0);
  }

  if(ileft < imid) { // a[01]
    if(size == 2 * _K_) { // leaf
      bv_pb(&(a->t), 1);
      bv_pb(&(a->t), 1);
      bv_pb(&(a->t), 0);
      bv_pb(&(a->t), 0);
      uint8_t leaf = encode_leaf(ia + ileft, imid - ileft, left);
      k2bp_write_leaf(a, leaf);
    } else {
      reccreate_k2bp(ia + ileft, imid - ileft, left, size / 2, a);
    }
  } else {
    bv_pb(&(a->t), 1);
    bv_pb(&(a->t), 0);
  }

  if(iright > imid) { // a[10]
    if(size == 2 * _K_) { // leaf
      bv_pb(&(a->t), 1);
      bv_pb(&(a->t), 1);
      bv_pb(&(a->t), 0);
      bv_pb(&(a->t), 0);
      uint8_t leaf = encode_leaf(ia + imid, iright - imid, mid);
      k2bp_write_leaf(a, leaf);
    } else {
      reccreate_k2bp(ia + imid, iright - imid, mid, size / 2, a);
    }
  } else {
    bv_pb(&(a->t), 1);
    bv_pb(&(a->t), 0);
  }

  if(iright < n) { // a[11]
    if(size == 2 * _K_) { // leaf
      bv_pb(&(a->t), 1);
      bv_pb(&(a->t), 1);
      bv_pb(&(a->t), 0);
      bv_pb(&(a->t), 0);
      uint8_t leaf = encode_leaf(ia + iright, n - iright, right);
      k2bp_write_leaf(a, leaf);
    } else {
      reccreate_k2bp(ia + iright, n - iright, right, size / 2, a);
    }
  } else {
    bv_pb(&(a->t), 1);
    bv_pb(&(a->t), 0);
  }

  bv_pb(&(a->t), 0); // finishing
}

static size_t create_k2bp(uint64_t ia[], size_t n, size_t rmsize, k2bp_t *a) {
  assert(ia != NULL && a != NULL);
  assert(n > 0);
  assert(rmsize > 1);
  k2bp_free(a);
  bv_init(&(a->t));
  a->l = (uint8_t*) malloc(sizeof(uint8_t) * 10);
  a->maxn_l = 10;
  a->n_l = 0;
  a->rmsize = rmsize;
  a->msize = k2bp_compute_height(rmsize);
  a->m = n;

  reccreate_k2bp(ia, n, 0, a->msize, a);

  a->maxn_l = 2 * ((a->n_l + 1) / 2);
  a->l = (uint8_t*) realloc(a->l, sizeof(uint8_t) * a->maxn_l / 2);
  bv_shrink(&(a->t));
  return a->t.n;
}

// the following functions are originally made by Giovanni Manzini
// source: https://github.com/acubeLab/k2tree/blob/main/k2text.c#L856

// interleaves two 32 bits integers in a single uint64_t 
// the bits of a (row index) are more significant than
// those of b (column index) because of how we number submatrices
static uint64_t bits_interleave(int64_t a, int64_t b) {
  uint64_t r = 0;
  assert(a<=UINT32_MAX && b <= UINT32_MAX);
  int c = 0;
  while(a!=0 || b!=0) {
    r |= (b&1)<<c++;
    r |= (a&1)<<c++;
    a >>= 1; b>>=1;  
    assert(c<=64);
  }
  return r;
}

static int uint64_cmp(const void *p, const void *q) {
  const uint64_t *a = p;
  const uint64_t *b = q;

  if(*a < *b) return -1;
  else if(*a > *b) return 1;
  return 0;
}

// create and return an interleaved array from the list of entries in a text file
// the matrix size stored in :msize is computed as follows: 
//  if xsize==0 *msize = largest index + 1
//  if xsize>0 that value is forced to be the matrix size (all indexes must be <xsize)
// since entries are encoded in 64 bits, each index can be at most 32 bits
// so the maximum matrix size is 2^32 (change ia[] type to go further)
static uint64_t *create_ia(FILE *f, size_t *n, size_t *msize, size_t xsize) {
  int64_t maxentry = 0; // largest entry in the file
  size_t size=10;      // current size of ia[]
  size_t i=0;          // elements in ia[]
  uint64_t *ia = malloc(size*sizeof(*ia));
  if(ia==NULL) quit("create_ia: malloc failed",__LINE__,__FILE__);

  int64_t a,b; size_t line=0;  
  while(1) {
    line++;
    int e = fscanf(f,"%" SCNd64 " %" SCNd64,&a,&b);
    if(e==EOF) break;
    // check input
    if(e!=2) {
      fprintf(stderr,"Invalid file content at line %zu\n",line);
      exit(EXIT_FAILURE);
    }
    if(a<0 || b<0) {
      fprintf(stderr,"Negative index at line %zu\n",line);
      exit(EXIT_FAILURE);
    }
    // since we are storing entries in 64 bits each index must fit in 32 bits       
    if(a>UINT32_MAX || b>UINT32_MAX) {
      fprintf(stderr,"Index too large at line %zu\n",line);
      exit(EXIT_FAILURE);
    }
    if(xsize>0 && (a>=xsize || b>=xsize)) {
      fprintf(stderr,"Index larger than the assigned size at line %zu\n",line);
      exit(EXIT_FAILURE);
    }
    // update maxentry
    if(a>maxentry) maxentry=a;
    if(b>maxentry) maxentry=b;
    // compute interleaved value
    uint64_t entry = bits_interleave(a,b);
    // enlarge ia if necessary
    if(i==size) {
      size = size*2;
      ia = realloc(ia,size*sizeof(*ia));
      if(ia==NULL) quit("create_ia: realloc failed",__LINE__,__FILE__);
    }
    assert(size>i);
    ia[i++] = entry;
  }
  // final resize
  size = i;
  ia = realloc(ia,size*sizeof(*ia));
  if(ia==NULL) quit("create_ia: realloc failed",__LINE__,__FILE__);
  // sort interleaved entries
  qsort(ia, size, sizeof(*ia), &uint64_cmp);
  // save output parameters   
  if(xsize==0) { // if xsize==0 size is largest index + 1
    if(maxentry+1>SIZE_MAX)  // highly unlikely, but you never know... 
      quit("create_ia: cannot represent matrix size",__LINE__,__FILE__);
    *msize = (size_t) maxentry+1;
  }
  else {  // if parameter xsize>0 that is the desired matrix size
    assert(maxentry<xsize);
    *msize = xsize;
  }
  *n = size;
  return ia;  
}

static size_t binsearch(uint64_t *ia, size_t n, uint64_t x) {
  assert(ia!=NULL && n>0);
  size_t l=0, r=n-1;
  while(l<r) {
    size_t m = (l+r)/2;
    if(ia[m]<x) l=m+1;
    else if(ia[m]==x) return m;
    else r=m; // replace with r = m-1 and later return r+1?
  }
  assert(l==r);
  if(ia[l]<x) {
    assert(r==n-1);
    return n;   // replace with return r+1?
  }
  return l;
}

// write error message and exit
static void quit(const char *msg, int line, char *file) {
  if(errno==0)  fprintf(stderr,"== %d == %s\n",getpid(), msg);
  else fprintf(stderr,"== %d == %s: %s\n",getpid(), msg,
      strerror(errno));
  fprintf(stderr,"== %d == Line: %d, File: %s\n",getpid(),line,file);
  exit(1);
}
