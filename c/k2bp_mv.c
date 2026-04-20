#include <getopt.h>
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#include "k2bp/k2bp.h"

void usage_and_exit(char* argv);

#define EXT ".k2bp"

int main(int argc, char* argv[]) {
  extern char *optarg;
  extern int optind, opterr, optopt;

  int c;

  while((c=getopt(argc, argv, "h")) != -1) {
    switch(c) {
      case 'h':
        usage_and_exit(argv[0]);
      case '?':
        fprintf(stderr, "Unkown option: %s", optarg);
        exit(1);
    }
  }

  optind -= 1;
  if(argc - optind != 3) usage_and_exit(argv[0]);
  argv += optind; argc -= optind;

  char fname[1000];
  strcpy(fname, argv[1]);
  k2bp_t a = K2BP_INITIALIZER;
  k2bp_load_from_file(&a, fname);
  double* vec_test = (double*) malloc(sizeof(double) * a.rmsize);
  FILE* f;
  f = fopen(argv[2], "r");
  size_t bytes_read = fread(vec_test, sizeof(double), a.rmsize, f);
  if(bytes_read != a.rmsize) {
    fprintf(stderr, "error reading vector file\n");
    exit(1);
  }
  fclose(f);
  double not_opt = 0;
  for(size_t i = 0; i < 100; i++) {
    double* ret = (double*) malloc(sizeof(double) * a.rmsize);
    k2bp_mv(&a, vec_test, ret);
    not_opt += ret[0];
    free(vec_test); vec_test = NULL;
    vec_test = ret;
  }

  printf("%f\n", not_opt);

  k2bp_free(&a);

  return 0;
}

void usage_and_exit(char* name) {
  fprintf(stderr, "Usage:\n\t%s [options] filename1 filename2 \n\n", name);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "\t-h      show this help message\n");    
  fprintf(stderr, "Sum compressed matrices in filename1 and filename2\n\n");
  exit(1);
}
