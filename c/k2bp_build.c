#include <getopt.h>
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#include "csr/csr.h"
#include "k2bp/k2bp.h"
#include "util/cmp.h"

void usage_and_exit(char* argv);

#define EXT ".k2bp"

int main(int argc, char* argv[]) {
  extern char *optarg;
  extern int optind, opterr, optopt;

  float p = -1LL;
  size_t size = 0;
  int check = 0;
  int cleaves = 0;
  int not_save = 0;
  int c;

  while((c=getopt(argc, argv, "t:cs:hn")) != -1) {
    switch(c) {
      case 'l':
        cleaves = 1; break;
      case 'c':
        check = 1; break;
      case 's':
        size = atoll(optarg); break;
      case 'n':
        not_save = 1;
      case 'h':
        usage_and_exit(argv[0]);
      case '?':
        fprintf(stderr, "Unkown option: %s", optarg);
        exit(1);
    }
  }

  optind -= 1;
  if(argc - optind != 2) usage_and_exit(argv[0]);
  argv += optind; argc -= optind;
  // checking values of args
  if(p != -1LL) {
    if(p < 0) {
      fprintf(stderr, "percentage of threshold (-t %f) should not be negative\n", p);
      exit(1);
    }
  }

  char fname[1000];
  strcpy(fname, argv[1]);
  k2bp_t a = K2BP_INITIALIZER;
  k2bp_build_from_textfile(&a, fname, size);
  if(cleaves) {
    k2bp_compress_leaves(&a);
  }

  if(not_save == 0) {
    char k2fname[1000];
    strcpy(k2fname, fname);
    strcat(k2fname, EXT);
    k2bp_save_to_file(&a, k2fname);
  }

  if(check) {
    csr_t matrix;
    csr_read_from_textfile(&matrix, argv[1], size);

    size_t n_matrix;
    uint32_t* ones_matrix = csr_nonzeros(&matrix, &n_matrix);

    qsort(ones_matrix, n_matrix, sizeof(uint32_t), &cmp);

    size_t n_k2bp;
    uint32_t* ones_k2bp = k2bp_nonzeros(&a, &n_k2bp);
    qsort(ones_k2bp, n_k2bp, sizeof(uint32_t), &cmp);

    if(n_matrix != n_k2bp) {
      fprintf(stderr, "ERROR COMPRESSION! amount of 1's doesn't match. expected: %zu got: %zu\n", n_matrix, n_k2bp);
      exit(1);
    }
    for(size_t i = 0; i < n_k2bp; i += 2) {
      if(ones_k2bp[i] != ones_matrix[i] || ones_k2bp[i + 1] != ones_matrix[i + 1]) {
        fprintf(stderr, "ERROR COMPRESSION! %zu-th 1 doesn't match. expected: (%" PRIu32 ", %" PRIu32 ") got: (%" PRIu32 ", %" PRIu32 ")\n", i, ones_matrix[i], ones_matrix[i + 1], ones_k2bp[i], ones_k2bp[i + 1]);
        exit(1);
      }
    }
    free(ones_matrix);
    free(ones_k2bp);
    csr_free(&matrix);
  }

  k2bp_free(&a);

  return 0;
}

void usage_and_exit(char* name) {
  fprintf(stderr, "Usage:\n\t%s [options] filename \n\n", name);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "\t-c      compressed->decompress->check\n");
  fprintf(stderr, "\t-l      compress leaves using rrr vector for better space\n");
  fprintf(stderr, "\t-s S    matrix actual size (def. largest index + 1)\n");
  fprintf(stderr, "\t-h      show this help message\n");    
  fprintf(stderr, "\t-n      don't save matrix, just check\n");    
  fprintf(stderr, "Compress filename\n\n");
  exit(1);
}
