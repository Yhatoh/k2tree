CC = gcc
CFLAGS  = -O3
LDFLAGS =
UTIL_DIR = util
K2BP_DIR = k2bp
CSR_DIR  = csr
LIBSAIS  = ../libsais/liblibsais.a

UTIL_SRC   = $(wildcard $(UTIL_DIR)/*.c)
K2BP_SRC   = $(wildcard $(K2BP_DIR)/*.c)
CSR_SRC    = $(wildcard $(CSR_DIR)/*.c)

debug: CFLAGS = -g -O0
debug: all

bv_test: bv.c $(UTIL_SRC)
	$(CC) $^ $(CFLAGS) -o bv.x

k2bp_build: k2bp_build.c $(UTIL_SRC) $(K2BP_SRC) $(CSR_SRC) $(LIBSAIS)
	$(CC) $^ $(CFLAGS) -o k2bp_build.x

k2bp_info: k2bp_info.c $(UTIL_SRC) $(K2BP_SRC) $(LIBSAIS)
	$(CC) $^ $(CFLAGS) -o k2bp_info.x

k2bp_addinfo: k2bp_addinfo.c $(UTIL_SRC) $(K2BP_SRC) $(LIBSAIS)
	$(CC) $^ $(CFLAGS) -lm -o k2bp_addinfo.x

k2bp_sum: k2bp_sum.c $(UTIL_SRC) $(K2BP_SRC) $(LIBSAIS)
	$(CC) $^ $(CFLAGS) -lm -o k2bp_sum.x

k2bp_mul: k2bp_mul.c $(UTIL_SRC) $(K2BP_SRC) $(LIBSAIS)
	$(CC) $^ $(CFLAGS) -lm -o k2bp_mul.x

k2bp_compress: k2bp_compress.c $(UTIL_SRC) $(K2BP_SRC) $(LIBSAIS)
	$(CC) $^ $(CFLAGS) -lm -o k2bp_compress.x

genrandmatrix: genrandmatrix.c
	$(CC) $^ $(CFLAGS) -o genrandmatrix.x

all: bv_test k2bp_build k2bp_info k2bp_addinfo k2bp_sum k2bp_mul genrandmatrix k2bp_compress

release: CFLAGS += -DNDEBUG
release: all

clean:
	rm -f *.x
