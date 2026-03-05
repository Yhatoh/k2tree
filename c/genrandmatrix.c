#include <assert.h>
#include <getopt.h>
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

void usage_and_exit(char* argv);

#define DEFAULTNAME "matrix.txt"

int main(int argc, char* argv[]) {
  extern char *optarg;
  extern int optind, opterr, optopt;

  uint32_t size = 100;
  float density = 0.2;
  int c;

  while((c=getopt(argc, argv, "hs:d:")) != -1) {
    switch(c) {
      case 's':
        size = atoi(optarg); break;
      case 'd':
        density = atof(optarg); break;
      case 'h':
        usage_and_exit(argv[0]);
      case '?':
        fprintf(stderr, "Unkown option: %s", optarg);
        exit(1);
    }
  }

  optind -= 1;
  if(!(argc - optind >= 1 && argc - optind <= 2)) usage_and_exit(argv[0]);
  argv += optind; argc -= optind;

  uint32_t m = size * size * density;
  uint32_t total = size * size;

  uint32_t *indices = malloc(total * sizeof(uint32_t));

  for (uint32_t i = 0; i < total; i++) {
    indices[i] = i;
  }

  fprintf(stdout, "creating file with name %s\n", argv[1]);

  FILE* f;
  if(argc == 1)
    f = fopen(DEFAULTNAME, "w");
  else
    f = fopen(argv[1], "w");


  uint32_t coords[m * 2];
  uint32_t written = 0;
  for (uint32_t i = 0; i < m; i++) {
    uint32_t j = i + rand() % (total - i);

    uint32_t temp = indices[i];
    indices[i] = indices[j];
    indices[j] = temp;

    uint32_t fila = indices[i] / size;
    uint32_t col = indices[i] % size;
    fprintf(f, "%" PRIu32 " %" PRIu32 "\n", fila, col);

    assert(written < 2 * m);
    coords[written] = fila;
    coords[written + 1] = col;
    written += 2;
  }

  fclose(f);

  for(size_t i = 0; i < m * 2; i += 2) {
    for(size_t j = 0; j < m * 2; j += 2) {
      if(i == j) continue;

      if(coords[i] == coords[j] &&
         coords[i + 1] == coords[j + 1]) {
        fprintf(stderr, "ERROR! Duplicate\n"); // this should never happen
        exit(1);
      }
    }
  }

  free(indices);
  return 0;
}

void usage_and_exit(char* name) {
  fprintf(stderr, "Usage:\n\t%s [options] (foutput) \n\n", name);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "\t-s S    size of the desired matrix (def. 100)\n");    
  fprintf(stderr, "\t-d D    density of the matrix (def. 0.2)\n");    
  fprintf(stderr, "\t-h      show this help message\n");    
  fprintf(stderr, "Generate a random SxS matrix of density D and is stored in foutput\n\n");
  exit(1);
}
