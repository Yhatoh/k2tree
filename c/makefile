bv_test:
	gcc bv.c k2bp/bv_t.c k2bp/randomer_t.c -O3 -o bv.x

k2bp_build:
	gcc k2bp_build.c k2bp/bv_t.c k2bp/k2bp.c csr/csr.c util/cmp.c -O3 -o k2bp_build.x

k2bp_info:
	gcc k2bp_info.c k2bp/bv_t.c k2bp/k2bp.c -O3 -o k2bp_info.x

k2bp_addinfo:
	gcc k2bp_addinfo.c k2bp/bv_t.c k2bp/k2bp.c -O3 -o k2bp_addinfo.x

all: bv_test k2bp_build k2bp_info k2bp_addinfo
