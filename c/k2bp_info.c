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
  if(argc - optind != 2) usage_and_exit(argv[0]);
  argv += optind; argc -= optind;

  char fname[1000];
  strcpy(fname, argv[1]);
  k2bp_t a = K2BP_INITIALIZER;
  k2bp_load_from_file(&a, fname);
  k2bp_show_stats(&a, fname, stdout);
  k2bp_free(&a);

  return 0;
}

void usage_and_exit(char* name) {
  fprintf(stderr, "Usage:\n\t%s [options] filename \n\n", name);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "\t-h      show this help message\n");    
  fprintf(stderr, "Show info about compress matrix in filename\n\n");
  exit(1);
}
