# Default directories for SDSL and LIBSAIS; can be overridden from the command line.
SDSL_DIR ?= ~/sdsl-lite
LIBSAIS_DIR ?= libsais

# Default flags; if "release" is a goal, set release flags.
FLAGS ?=
ifneq ($(findstring release,$(MAKECMDGOALS)),)
FLAGS := -O3 -m64 -DNDEBUG
endif

k2bp_build:
	g++ -I $(SDSL_DIR)/include -L $(SDSL_DIR)/lib k2bp_build.cpp $(LIBSAIS_DIR)/liblibsais.a -lsdsl -ldivsufsort -ldivsufsort64 $(FLAGS) -o k2bp_build.x

k2bp_comp_build:
	g++ -I $(SDSL_DIR)/include -L $(SDSL_DIR)/lib k2bp_compr.cpp $(LIBSAIS_DIR)/liblibsais.a -lsdsl -ldivsufsort -ldivsufsort64 $(FLAGS) -o k2bp_compr.x

k2bp_info:
	g++ -I $(SDSL_DIR)/include -L $(SDSL_DIR)/lib k2bp_info.cpp $(LIBSAIS_DIR)/liblibsais.a -lsdsl -ldivsufsort -ldivsufsort64 $(FLAGS) -DINFO_SPACE -o k2bp_info.x

k2bp_comp_info:
	g++ -I $(SDSL_DIR)/include -L $(SDSL_DIR)/lib k2bp_compr_info.cpp $(LIBSAIS_DIR)/liblibsais.a -lsdsl -ldivsufsort -ldivsufsort64 $(FLAGS) -DINFO_SPACE -o k2bp_compr_info.x

k2bp_mult_1:
	g++ -I $(SDSL_DIR)/include -L $(SDSL_DIR)/lib k2bp_mult_1.cpp $(LIBSAIS_DIR)/liblibsais.a -lsdsl -ldivsufsort -ldivsufsort64 $(FLAGS) -o k2bp_mult_1.x

k2bp_mult_comp_1:
	g++ -I $(SDSL_DIR)/include -L $(SDSL_DIR)/lib k2bp_compr_mult.cpp $(LIBSAIS_DIR)/liblibsais.a -lsdsl -ldivsufsort -ldivsufsort64 $(FLAGS) -o k2bp_compr_mult.x

k2bp_rand_test:
	g++ -I $(SDSL_DIR)/include -L $(SDSL_DIR)/lib test.cpp $(LIBSAIS_DIR)/liblibsais.a -lsdsl -ldivsufsort -ldivsufsort64 $(FLAGS) -o test.x

gen_rand_matrix:
	g++ gen_rand_matrix.cpp $(FLAGS) -o gen_rand_matrix.x

all: k2bp_build k2bp_comp_build k2bp_info k2bp_comp_info k2bp_mult_1 k2bp_mult_comp_1 k2bp_rand_test gen_rand_matrix

release: all

