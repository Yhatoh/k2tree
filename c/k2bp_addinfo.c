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

  int subinfo = 0;
  int excinfo = 0;
  float p = 1;
  int c;

  while((c=getopt(argc, argv, "hesp:")) != -1) {
    switch(c) {
      case 'e':
        excinfo = 1; break;
      case 's':
        subinfo = 1; break;
      case 'p':
        p = atof(optarg); break;
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

  if(subinfo) {
    k2bp_addsubtree_info(&a, sqrt(a.t.n / 2) * p);
    k2bp_checksubtree_info(&a);
  }

  if(excinfo) {
    k2bp_build_exc_sampling(&a);
  }

  k2bp_save_to_file(&a, fname);

  k2bp_free(&a);

  return 0;
}

void usage_and_exit(char* name) {
  fprintf(stderr, "Usage:\n\t%s [options] filename \n\n", name);
  fprintf(stderr, "Options:\n");

  fprintf(stderr, "\t-e      add sampling of excess\n");    
  fprintf(stderr, "\t-s      add subtree information\n");    
  fprintf(stderr, "\t-p P    only save subtree information with at least P * sqrt(N) nodes (def. P = 1)\n");    
  fprintf(stderr, "\t-h      show this help message\n");    
  fprintf(stderr, "Add information for speedup to compressed matrix in filename\n\n");
  exit(1);
}
