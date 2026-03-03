#include <assert.h>
#include <errno.h>
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "csr.h"

static void quit(const char *msg, int line, char *file);

void csr_read_from_textfile(csr_t *a, char *fname, size_t fsize) {
  assert(a != NULL && fname != NULL);
  FILE* f = fopen(fname, "rt");
  if(f == NULL) quit("csr_read_from_textfile: file cannot be open", __LINE__, __FILE__);

  uint32_t n = 0;
  uint32_t maxn = 10;
  uint32_t* ones = (uint32_t*) malloc(sizeof(uint32_t) * maxn);
  int64_t x, y;
  size_t maxindex;
  size_t lines = 0;
  while(1) {
    lines++;
    int e = fscanf(f,"%" SCNd64 " %" SCNd64,&x,&y); 
    if(e == EOF)
      break;
    if(e != 2) {
      char message[1000];
      sprintf(message, "csr_read_from_textfile: invalid content at line %zu", lines);
      quit(message, __LINE__, __FILE__);
    }
    if(x < 0 || y < 0) {
      char message[1000];
      sprintf(message, "csr_read_from_textfile: negative index at line %zu", lines);
      quit(message, __LINE__, __FILE__);
    }
    if(x > UINT32_MAX || y > UINT32_MAX) {
      char message[1000];
      sprintf(message, "csr_read_from_textfile: index at line %zu bigger than assigned size", lines);
      quit(message, __LINE__, __FILE__);
    }

    if(x > maxindex) maxindex = x;
    if(y > maxindex) maxindex = y;
    if(n >= maxn) {
      maxn *= 2;
      ones = (uint32_t*) realloc(ones, sizeof(uint32_t) * maxn);
      if(ones == NULL)
        quit("csr_read_from_textfile: realloc failed", __LINE__, __FILE__);
    }

    ones[n] = x;
    ones[n + 1] = y;
    n += 2;
  }
  fclose(f);

  maxn = n;
  ones = (uint32_t*) realloc(ones, sizeof(uint32_t) * maxn);
  if(ones == NULL)
    quit("csr_read_from_textfile: realloc failed", __LINE__, __FILE__);

  a->n = n;
  a->nonzeros = ones;
  if(fsize == 0) {
    a->msize = maxindex + 1;
  } else {
    a->msize = fsize;
  }
}

void csr_free(csr_t* a) {
  if(a->nonzeros != NULL) {
    free(a->nonzeros);
    a->nonzeros = NULL;
  }
  a->n = 0;
  a->msize = 0;
}

uint32_t* csr_nonzeros(csr_t* a, size_t* n) {
  uint32_t* arr = (uint32_t*) malloc(sizeof(uint32_t) * a->n);
  for(size_t i = 0; i < a->n; i += 2) {
    arr[i] = a->nonzeros[i];
    arr[i + 1] = a->nonzeros[i + 1];
  }
  *n = a->n;
  return arr;
}

// write error message and exit
static void quit(const char *msg, int line, char *file) {
  if(errno==0)  fprintf(stderr,"== %d == %s\n",getpid(), msg);
  else fprintf(stderr,"== %d == %s: %s\n",getpid(), msg,
               strerror(errno));
  fprintf(stderr,"== %d == Line: %d, File: %s\n",getpid(),line,file);
  exit(1);
}
