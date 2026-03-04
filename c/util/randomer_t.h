#ifndef __RANDOMER_H__
#define __RANDOMER_H__

#include <stdint.h>
#include <stddef.h>
#include <time.h>

typedef struct randomer_t {
  size_t min, max;
  size_t seed;
} randomer_t;

void randomer_t_init(randomer_t* r, size_t min, size_t max, size_t seed);
size_t randomer_t_gennum(randomer_t* r);

#endif
