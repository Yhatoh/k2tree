#include "k2bp/bv_t.h"
#include "k2bp/randomer_t.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>


uint8_t test_pb() {
  static const size_t N = 100000;
  uint8_t exp[N];
  randomer_t r;
  randomer_t_init(&r, 0, 1, 42);

  for(size_t i = 0; i < N; i++)
    exp[i] = randomer_t_gennum(&r);

  bv_t z;
  bv_init(&z);

  for(size_t i = 0; i < N; i++)
    bv_pb(&z, exp[i]);

  if(z.n != N) {
    fprintf(stderr, "== %d == ERROR: size %zu bv_t doesn't match with expected %zu\n", getpid(), z.n, N);
    fprintf(stderr, "== %d == Line: %d, File: %s\n",getpid(),__LINE__,__FILE__);
    exit(1);
  }
  for(size_t i = 0; i < N; i++) {
    if(bv_i(&z, i) != exp[i]) {
      fprintf(stderr, "== %d == ERROR: %zu-th element doesn't match: expected: %zu, got: %zu\n", getpid(), i, (size_t) exp[i], bv_i(&z, i));
      fprintf(stderr, "== %d == Line: %d, File: %s\n",getpid(),__LINE__,__FILE__);
      exit(1);
    }
  }
  bv_free(&z);
  return 1;
}

uint8_t test_append_int() {
  static const size_t N = 100;
  static const size_t Ni = 10000;
  uint64_t exp[N + Ni];
  randomer_t r;
  randomer_t_init(&r, 0, 1, 42);

  for(size_t i = 0; i < N; i++)
    exp[i] = randomer_t_gennum(&r);

  bv_t z;
  bv_init(&z);

  for(size_t i = 0; i < N; i++)
    bv_pb(&z, exp[i]);

  if(z.n != N) {
    fprintf(stderr, "== %d == ERROR: size %zu bv_t doesn't match with expected %zu\n", getpid(), z.n, N);
    fprintf(stderr, "== %d == Line: %d, File: %s\n",getpid(),__LINE__,__FILE__);
    exit(1);
  }
  for(size_t i = 0; i < N; i++) {
    if(bv_i(&z, i) != exp[i]) {
      fprintf(stderr, "== %d == ERROR: %zu-th element doesn't match: expected: %zu, got: %zu\n", getpid(), i, (size_t) exp[i], bv_i(&z, i));
      fprintf(stderr, "== %d == Line: %d, File: %s\n",getpid(),__LINE__,__FILE__);
      exit(1);
    }
  }

  randomer_t r2;
  randomer_t_init(&r2, 0, UINT64_MAX - 1, 42);
  for(size_t i = N; i < N + Ni; i++) {
    exp[i] = randomer_t_gennum(&r2);
    bv_append_int(&z, exp[i]);
  }

  if(z.n != N + Ni * 64) {
    fprintf(stderr, "== %d == ERROR: size %zu bv_t doesn't match with expected %zu\n", getpid(), z.n, N + Ni * 64);
    fprintf(stderr, "== %d == Line: %d, File: %s\n",getpid(),__LINE__,__FILE__);
    exit(1);
  }
  for(size_t i = 0; i < Ni; i++) {
    if(bv_get_int(&z, (i * 64) + N, 64) != exp[N + i]) {
      fprintf(stderr, "== %d == ERROR: %zu-th element doesn't match: expected: %zu, got: %zu\n", getpid(), i * 64 + N, (size_t) exp[i + N], bv_get_int(&z, (i * 64) + N, 64));
      fprintf(stderr, "== %d == Line: %d, File: %s\n",getpid(),__LINE__,__FILE__);
      exit(1);
    }
  }

  bv_free(&z);
  return 1;
}

uint8_t test_append() {
  static const size_t N = 10000;
  static const size_t Ni = 10000;
  uint64_t exp[N + Ni];
  randomer_t r;
  randomer_t_init(&r, 0, 1, 42);

  for(size_t i = 0; i < N + Ni; i++)
    exp[i] = randomer_t_gennum(&r);

  bv_t z;
  bv_init(&z);

  for(size_t i = 0; i < N; i++)
    bv_pb(&z, exp[i]);

  if(z.n != N) {
    fprintf(stderr, "== %d == ERROR: size %zu bv_t doesn't match with expected %zu\n", getpid(), z.n, N);
    fprintf(stderr, "== %d == Line: %d, File: %s\n",getpid(),__LINE__,__FILE__);
    exit(1);
  }
  for(size_t i = 0; i < N; i++) {
    if(bv_i(&z, i) != exp[i]) {
      fprintf(stderr, "== %d == ERROR: %zu-th element doesn't match: expected: %zu, got: %zu\n", getpid(), i, (size_t) exp[i], bv_i(&z, i));
      fprintf(stderr, "== %d == Line: %d, File: %s\n",getpid(),__LINE__,__FILE__);
      exit(1);
    }
  }
  bv_t z2;
  bv_init(&z2);

  for(size_t i = N; i < N + Ni; i++)
    bv_pb(&z2, exp[i]);

  if(z2.n != Ni) {
    fprintf(stderr, "== %d == ERROR: size %zu bv_t doesn't match with expected %zu\n", getpid(), z2.n, Ni);
    fprintf(stderr, "== %d == Line: %d, File: %s\n",getpid(),__LINE__,__FILE__);
    exit(1);
  }
  for(size_t i = N; i < Ni; i++) {
    if(bv_i(&z2, i - N) != exp[i]) {
      fprintf(stderr, "== %d == ERROR: %zu-th element doesn't match: expected: %zu, got: %zu\n", getpid(), i, (size_t) exp[i], bv_i(&z2, i - N));
      fprintf(stderr, "== %d == Line: %d, File: %s\n",getpid(),__LINE__,__FILE__);
      exit(1);
    }
  } 

  bv_append(&z, &z2);
  for(size_t i = 0; i < N + Ni; i++) {
    if(bv_i(&z, i) != exp[i]) {
      fprintf(stderr, "== %d == ERROR: %zu-th element doesn't match: expected: %zu, got: %zu\n", getpid(), i, (size_t) exp[i], bv_i(&z, i));
      fprintf(stderr, "== %d == Line: %d, File: %s\n",getpid(),__LINE__,__FILE__);
      exit(1);
    }
  }

  bv_free(&z);
  bv_free(&z2);
  return 1;
}

// in the future add arguments and bla bla
int main() {
  srand(42);
  test_pb();
  test_append_int();
  test_append();
  return 0;
}
