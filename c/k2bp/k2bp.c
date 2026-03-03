#include <assert.h>
#include <errno.h>
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "k2bp.h"

static size_t binsearch(uint64_t *ia, size_t n, uint64_t x);
static uint64_t *create_ia(FILE *f, size_t *n, size_t *msize, size_t xsize);
static size_t create_k2bp(uint64_t ia[], size_t n, size_t msize, k2bp_t *a);
static void quit(const char *msg, int line, char *file);
static uint8_t encode_leaf(uint64_t ia[], size_t n, size_t smin);
static void reck2bp_nonzeros(k2bp_traversal_t* pos_a, const k2bp_t* a, uint32_t* arr, size_t* n);

size_t k2bp_build_from_textfile(k2bp_t *a, char* fname, size_t xsize) {
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

uint8_t k2bp_read_leaf(const k2bp_t*a, size_t pos) {
  assert(pos < a->n_l);
  if(pos % 2 == 0)
    return a->l[pos / 2] & 15;
  return (a->l[pos / 2] >> 4) & 15;
}

size_t k2bp_compute_height(size_t rmsize) {
  assert(rmsize > 1);
  size_t msize = 2;
  // task for future me:
  //  add detection of overflow
  while(msize < rmsize) msize *= 2;
  return msize;
}

void k2bp_free(k2bp_t* a) {
  a->msize = a->rmsize = a->m = a->maxn_l = a->n_l = 0;
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

  a->l = NULL;
  a->exc_min_samples = NULL;
  a->exc_samples = NULL;
  a->leaves_samples = NULL;
  a->subtreeinfo = NULL;
  a->leavesinfo = NULL;
}

uint32_t* k2bp_nonzeros(const k2bp_t* a, size_t* n) {
  uint32_t* arr = (uint32_t*) malloc(sizeof(uint32_t) * (a->m * 2));
  *n = 0;
  k2bp_traversal_t pos_a = {a->msize, 0, 0, 0, 0};
  reck2bp_nonzeros(&pos_a, a, arr, n);
  return arr;
}

void k2bp_save_to_file(const k2bp_t* a, const char* fname) {
  assert(a != NULL);
  assert(a->l != NULL && a->t.a != NULL);

  char info_name[1000];
  strcpy(info_name, fname);
  strcat(info_name, ".i");

  FILE* f = fopen(info_name, "w");
  if(f == NULL)
    quit("k2bp_write_to_file: file cannot be open", __LINE__, __FILE__);

  size_t w = fwrite(&(a->msize), sizeof(size_t), 1, f);
  printf("%zu %zu\n", w, sizeof(size_t));
  if(w != 1)
    quit("k2bp_write_to_file: error writing in file", __LINE__, __FILE__);
  w = fwrite(&(a->rmsize), sizeof(size_t), 1, f);
  if(w != 1)
    quit("k2bp_write_to_file: error writing in file", __LINE__, __FILE__);
  w = fwrite(&(a->m), sizeof(size_t), 1, f);
  if(w != 1)
    quit("k2bp_write_to_file: error writing in file", __LINE__, __FILE__);

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
    quit("k2bp_write_to_file: file cannot be open", __LINE__, __FILE__);

  w = fwrite(&(a->maxn_l), sizeof(size_t), 1, f);
  if(w != 1)
    quit("k2bp_write_to_file: error writing in file", __LINE__, __FILE__);
  w = fwrite(&(a->n_l), sizeof(size_t), 1, f);
  if(w != 1)
    quit("k2bp_write_to_file: error writing in file", __LINE__, __FILE__);
  w = fwrite(a->l, sizeof(uint8_t), a->n_l, f);
  if(w != a->n_l)
    quit("k2bp_write_to_file: error writing in file", __LINE__, __FILE__);

  fclose(f);

  if(a->exc_min_samples != NULL) {
    char exc_name[1000];
    strcpy(exc_name, fname);
    strcat(exc_name, ".exc");

    f = fopen(exc_name, "w");
    if(f == NULL)
      quit("k2bp_write_to_file: file cannot be open", __LINE__, __FILE__);

    w = fwrite(&(a->exc_min_samples), sizeof(uint16_t), (a->t.n + BLOCK_SIZE - 1) / BLOCK_SIZE + 1, f);
    if(w != (a->t.n + BLOCK_SIZE - 1) / BLOCK_SIZE + 1)
      quit("k2bp_write_to_file: error writing in file", __LINE__, __FILE__);
    w = fwrite(&(a->exc_samples), sizeof(uint16_t), (a->t.n + BLOCK_SIZE - 1) / BLOCK_SIZE + 1, f);
    if(w != (a->t.n + BLOCK_SIZE - 1) / BLOCK_SIZE + 1)
      quit("k2bp_write_to_file: error writing in file", __LINE__, __FILE__);
    w = fwrite(&(a->leaves_samples), sizeof(uint8_t), (a->t.n + BLOCK_SIZE - 1) / BLOCK_SIZE + 1, f);
    if(w != (a->t.n + BLOCK_SIZE - 1) / BLOCK_SIZE + 1)
      quit("k2bp_write_to_file: error writing in file", __LINE__, __FILE__);

    fclose(f);
  }

  if(a->subtreeinfo != NULL) {
    char sub_name[100];
    strcpy(sub_name, fname);
    strcat(sub_name, ".s");

    f = fopen(sub_name, "w");
    if(f == NULL)
      quit("k2bp_write_to_file: file cannot be open", __LINE__, __FILE__);

    w = fwrite(&(a->n_info), sizeof(size_t), 1, f);
    if(w != 1)
      quit("k2bp_write_to_file: error writing in file", __LINE__, __FILE__);
    w = fwrite(&(a->subtreeinfo), sizeof(uint64_t), a->n_info, f);
    if(w != a->n_info)
      quit("k2bp_write_to_file: error writing in file", __LINE__, __FILE__);
    w = fwrite(&(a->leavesinfo), sizeof(uint64_t), a->n_info, f);
    if(w != a->n_info)
      quit("k2bp_write_to_file: error writing in file", __LINE__, __FILE__);

    fclose(f);
  }
}

void k2bp_load_from_file(k2bp_t* a, const char* fname) {
  assert(a != NULL);

  char info_name[1000];
  strcpy(info_name, fname);
  strcat(info_name, ".i");

  FILE* f = fopen(info_name, "r");
  if(f == NULL)
    quit("k2bp_write_to_file: file cannot be open", __LINE__, __FILE__);

  size_t w = fread(&(a->msize), sizeof(size_t), 1, f);
  if(w != 1)
    quit("k2bp_write_to_file: error writing in file", __LINE__, __FILE__);
  w = fread(&(a->rmsize), sizeof(size_t), 1, f);
  if(w != 1)
    quit("k2bp_write_to_file: error writing in file", __LINE__, __FILE__);
  w = fread(&(a->m), sizeof(size_t), 1, f);
  if(w != 1)
    quit("k2bp_write_to_file: error writing in file", __LINE__, __FILE__);

  fclose(f);

  char tree_name[1000];
  strcpy(tree_name, fname);
  strcat(tree_name, ".t");
  bv_save_to_file(&(a->t), tree_name);

  char l_name[1000];
  strcpy(l_name, fname);
  strcat(l_name, ".l");

  f = fopen(l_name, "r");
  if(f == NULL)
    quit("k2bp_write_to_file: file cannot be open", __LINE__, __FILE__);

  w = fread(&(a->maxn_l), sizeof(size_t), 1, f);
  if(w != 1)
    quit("k2bp_write_to_file: error writing in file", __LINE__, __FILE__);
  w = fread(&(a->n_l), sizeof(size_t), 1, f);
  if(w != 1)
    quit("k2bp_write_to_file: error writing in file", __LINE__, __FILE__);

  a->l = (uint8_t*) malloc(sizeof(uint8_t) * a->n_l);

  w = fread(a->l, sizeof(uint8_t), a->n_l, f);
  if(w != a->n_l)
    quit("k2bp_write_to_file: error writing in file", __LINE__, __FILE__);

  fclose(f);

  if(a->exc_min_samples != NULL) {
    char exc_name[1000];
    strcpy(exc_name, fname);
    strcat(exc_name, ".exc");

    f = fopen(exc_name, "r");
    if(f == NULL)
      quit("k2bp_write_to_file: file cannot be open", __LINE__, __FILE__);

    a->exc_min_samples = (uint16_t*) malloc(sizeof(uint16_t) * ((a->t.n + BLOCK_SIZE - 1) / BLOCK_SIZE + 1));
    w = fread(&(a->exc_min_samples), sizeof(uint16_t), (a->t.n + BLOCK_SIZE - 1) / BLOCK_SIZE + 1, f);
    if(w != (a->t.n + BLOCK_SIZE - 1) / BLOCK_SIZE + 1)
      quit("k2bp_write_to_file: error writing in file", __LINE__, __FILE__);

    a->exc_samples = (uint16_t*) malloc(sizeof(uint16_t) * ((a->t.n + BLOCK_SIZE - 1) / BLOCK_SIZE + 1));
    w = fread(&(a->exc_samples), sizeof(uint16_t), (a->t.n + BLOCK_SIZE - 1) / BLOCK_SIZE + 1, f);
    if(w != (a->t.n + BLOCK_SIZE - 1) / BLOCK_SIZE + 1)
      quit("k2bp_write_to_file: error writing in file", __LINE__, __FILE__);

    a->leaves_samples = (uint8_t*) malloc(sizeof(uint8_t) * ((a->t.n + BLOCK_SIZE - 1) / BLOCK_SIZE + 1));
    w = fread(&(a->leaves_samples), sizeof(uint8_t), (a->t.n + BLOCK_SIZE - 1) / BLOCK_SIZE + 1, f);
    if(w != (a->t.n + BLOCK_SIZE - 1) / BLOCK_SIZE + 1)
      quit("k2bp_write_to_file: error writing in file", __LINE__, __FILE__);

    fclose(f);
  }

  if(a->subtreeinfo != NULL) {
    char sub_name[100];
    strcpy(sub_name, fname);
    strcat(sub_name, ".s");

    f = fopen(sub_name, "r");
    if(f == NULL)
      quit("k2bp_write_to_file: file cannot be open", __LINE__, __FILE__);

    w = fread(&(a->n_info), sizeof(size_t), 1, f);
    if(w != 1)
      quit("k2bp_write_to_file: error writing in file", __LINE__, __FILE__);

    a->subtreeinfo = (uint64_t*) malloc(sizeof(uint64_t) * a->n_info);
    w = fread(&(a->subtreeinfo), sizeof(uint64_t), a->n_info, f);
    if(w != a->n_info)
      quit("k2bp_write_to_file: error writing in file", __LINE__, __FILE__);

    a->leavesinfo = (uint32_t*) malloc(sizeof(uint32_t) * a->n_info);
    w = fread(&(a->leavesinfo), sizeof(uint64_t), a->n_info, f);
    if(w != a->n_info)
      quit("k2bp_write_to_file: error writing in file", __LINE__, __FILE__);

    fclose(f);
  }
}

// ----------------------------------------------------------

// auxiliary functions
static void reck2bp_nonzeros(k2bp_traversal_t* pos_a, const k2bp_t* a, uint32_t* arr, size_t* n) {
  assert(pos_a->i_t < a->t.n && pos_a->i_l < a->n_l);
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

  k2bp_traversal_t pos_aux = {pos_a->msize / 2, pos_a->x, pos_a->y, pos_a->i_t + 1, pos_a->i_l};
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

  reccreate_k2bp(ia, n, 0, rmsize, a);

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
