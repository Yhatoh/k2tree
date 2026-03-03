#include <stdint.h>
#include "cmp.h"

int cmp(const void* a, const void* b) {
  uint32_t int_a = *((uint32_t*) a);
  uint32_t int_b = *((uint32_t*) b);

  if(int_a == int_b) return 0;
  else if(int_a < int_b) return -1;
  else return 1;
}
