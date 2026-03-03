#include <inttypes.h>
#include <stdlib.h>
#include "k2bp/k2bp.h"

int main(int argc, char* argv[]) {
  k2bp_t a = K2BP_INITIALIZER;
  k2bp_build_from_textfile(&a, argv[1], 16);

  size_t n;
  uint32_t* check = k2bp_nonzeros(&a, &n);
  uint32_t* expected = (uint32_t*) malloc(sizeof(uint32_t) * n);

  printf("%zu\n", a.m);
  int64_t lines = a.m;
  FILE* f = fopen(argv[1], "rt");
  for(size_t i = 0; i < lines * 2;i+=2) {
    int64_t x, y;
    fscanf(f,"%" SCNd64 " %" SCNd64,&x,&y); 
    expected[i] = x; expected[i + 1] = y;
  }
  fclose(f);

  for(size_t i = 0; i < n; i += 2) {
    printf("(%" PRIu32 ", %" PRIu32 ")\n", check[i], check[i + 1]);
  }
  
  printf("------------\n");
  for(size_t i = 0; i < n; i += 2) {
    printf("(%" PRIu32 ", %" PRIu32 ")\n", expected[i], expected[i + 1]);
  }
  k2bp_free(&a);
  free(expected);
  free(check);
  return 0;
}
