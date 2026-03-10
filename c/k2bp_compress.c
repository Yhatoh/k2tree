#include <assert.h>
#include <getopt.h>
#include <inttypes.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "k2bp/k2bp.h"

void usage_and_exit(char* argv);

#define EXT ".k2bp"

int main(int argc, char* argv[]) {
  extern char *optarg;
  extern int optind, opterr, optopt;

  int p = 32;
  int check = 0;
  int c;

  while((c=getopt(argc, argv, "hcp:")) != -1) {
    switch(c) {
      case 'c':
        check = 1; break;
      case 'p':
        p = atoi(optarg); break;
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

  char fname[1000];
  strcpy(fname, argv[1]);
  k2bp_t a = K2BP_INITIALIZER;
  k2bp_load_from_file(&a, fname);
  k2bp_t ca = K2BP_INITIALIZER;
  k2bp_compress_subtrees(&a, &ca, p);
  if(check) {
    k2bp_t da = K2BP_INITIALIZER;
    k2bp_decompress_subtrees(&ca, &da);
    assert(k2bp_equal(&da, &ca) == 1);
    k2bp_free(&da);
  }

  char fname_save[1000];
  strcpy(fname_save, argv[1]);
  strcat(fname_save, ".c");
  k2bp_save_to_file(&ca, fname_save);

  k2bp_free(&ca);
  k2bp_free(&a);

  return 0;
}

void usage_and_exit(char* name) {
  fprintf(stderr, "Usage:\n\t%s [options] filename \n\n", name);
  fprintf(stderr, "Options:\n");

  fprintf(stderr, "\t-p P    only compress subtrees with size at least P nodes (def. P = 32)\n");    
  fprintf(stderr, "\t-h      show this help message\n");    
  fprintf(stderr, "Compressed subtrees\n\n");
  exit(1);
}
