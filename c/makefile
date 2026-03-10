CC = gcc
CFLAGS = -O3
LDFLAGS =

bv_test:
	$(CC) bv.c util/bv_t.c util/randomer_t.c $(CFLAGS) -o bv.x

k2bp_build:
	$(CC) k2bp_build.c util/bv_t.c util/dsu.c util/vu64.c k2bp/k2bp.c csr/csr.c util/cmp.c $(CFLAGS) -o k2bp_build.x

k2bp_info:
	$(CC) k2bp_info.c util/bv_t.c util/dsu.c util/vu64.c k2bp/k2bp.c $(CFLAGS) -o k2bp_info.x

k2bp_addinfo:
	$(CC) k2bp_addinfo.c util/bv_t.c util/dsu.c util/vu64.c k2bp/k2bp.c $(CFLAGS) -lm -o k2bp_addinfo.x

k2bp_sum:
	$(CC) k2bp_sum.c util/bv_t.c util/dsu.c util/vu64.c k2bp/k2bp.c $(CFLAGS) -lm -o k2bp_sum.x

k2bp_mul:
	$(CC) k2bp_mul.c util/bv_t.c util/dsu.c util/vu64.c k2bp/k2bp.c $(CFLAGS) -lm -o k2bp_mul.x

k2bp_compress:
	$(CC) k2bp_compress.c util/bv_t.c util/dsu.c util/vu64.c k2bp/k2bp.c $(CFLAGS) -lm -o k2bp_compress.x

genrandmatrix:
	$(CC) genrandmatrix.c $(CFLAGS) -o genrandmatrix.x


all: bv_test k2bp_build k2bp_info k2bp_addinfo k2bp_sum k2bp_mul genrandmatrix k2bp_compress

release: CFLAGS += -DNDEBUG
release: all
