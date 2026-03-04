#ifndef _VU64DEFS_H
#define _VU64DEFS_H 1

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

// resizable vector of uint64's
// (at the moment the vector can only grow in size) 
typedef struct {
  uint64_t *v;       // array of u64's
  size_t n;          // number of elements 
  size_t nmax;       // maximum capacity 
} vu64_t;

// create small vector
void vu64_init(vu64_t *z);
// deallocate 
void vu64_free(vu64_t *z);
// add i elements at the end of z
void vu64_grow(vu64_t *z, size_t i);
// write to file
void vu64_write(FILE *f, vu64_t *z);

#endif /* _VU64DEFS_C */
