#ifndef __CSR_H_
#define __CSR_H_

#include <stddef.h>
#include <stdint.h>

typedef struct csr_t {
  size_t n; // amount of entries
  size_t msize;
  uint32_t* nonzeros;
} csr_t;

void csr_read_from_textfile(csr_t* a, char* fname, size_t fsize);
void csr_free(csr_t* a);
uint32_t* csr_nonzeros(csr_t* a, size_t* n);

#endif
