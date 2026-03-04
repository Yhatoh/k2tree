#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <stdbool.h>
#include <string.h>
#include <errno.h>
#include "vu64.h" 

static void quit(const char *msg, int line, char *file);

void vu64_init(vu64_t *z)
{
  z->nmax=16;
  z->n=0;
  z->v = malloc(z->nmax*sizeof *(z->v) );
}

void vu64_free(vu64_t *z)
{
  free(z->v);
  z->v = NULL; 
  z->nmax = z->n=0;
}

void vu64_write(FILE *f, vu64_t *z)
{
  size_t w = fwrite(z->v,sizeof(*(z->v)),z->n,f);
  if(w!=z->n) quit("Error writing vu64 to file",__LINE__,__FILE__);
}

// make sure there is space for i more elements at the end of z
void vu64_grow(vu64_t *z, size_t i) 
{
  z->n +=i;
  while(z->n>z->nmax) {
    z->nmax *= 2;
    z->v = realloc(z->v, z->nmax*sizeof *(z->v) );
    if(z->v==NULL) quit("realloc failed",__LINE__,__FILE__);
  }
}


// write error message and exit
static void quit(const char *msg, int line, char *file) {
  if(errno==0)  fprintf(stderr,"== %d == %s\n",getpid(), msg);
  else fprintf(stderr,"== %d == %s: %s\n",getpid(), msg,
               strerror(errno));
  fprintf(stderr,"== %d == Line: %d, File: %s\n",getpid(),line,file);
  exit(1);
}
