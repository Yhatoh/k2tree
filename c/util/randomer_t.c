#include <stdlib.h>
#include "randomer_t.h"

void randomer_t_init(randomer_t* r, size_t min, size_t max, size_t seed) {
  r->min = min;
  r->max = max;
  r->seed = seed;
}

size_t randomer_t_gennum(randomer_t* r) {
  return rand() % (r->max - r->min + 1) + r->min;
}
