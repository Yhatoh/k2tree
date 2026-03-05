bv_test:
	gcc bv.c util/bv_t.c util/randomer_t.c -O3 -o bv.x

k2bp_build:
	gcc k2bp_build.c util/bv_t.c util/vu64.c k2bp/k2bp.c csr/csr.c util/cmp.c -O3 -o k2bp_build.x

k2bp_info:
	gcc k2bp_info.c util/bv_t.c util/vu64.c k2bp/k2bp.c -O3 -o k2bp_info.x

k2bp_addinfo:
	gcc k2bp_addinfo.c util/bv_t.c util/vu64.c k2bp/k2bp.c -O3 -lm -o k2bp_addinfo.x

genrandmatrix:
	gcc genrandmatrix.c -O3 -o genrandmatrix.x


all: bv_test k2bp_build k2bp_info k2bp_addinfo genrandmatrix
