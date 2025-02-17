compile_test:
	g++ -I ~/local-sdsl-lite/include -L ~/local-sdsl-lite/lib test_file.cpp libsais/liblibsais.a -lsdsl -ldivsufsort -ldivsufsort64 -o out

k2bp_build:
	g++ -I ~/local-sdsl-lite/include -L ~/local-sdsl-lite/lib k2bp_build.cpp libsais/liblibsais.a -lsdsl -ldivsufsort -ldivsufsort64 -DNDEBUG -o k2bp_build.x;

k2bp_comp_build:
	g++ -I ~/local-sdsl-lite/include -L ~/local-sdsl-lite/lib k2bp_compr.cpp libsais/liblibsais.a -lsdsl -ldivsufsort -ldivsufsort64 -o k2bp_compr.x;

k2bp_info:
	g++ -I ~/local-sdsl-lite/include -L ~/local-sdsl-lite/lib k2bp_info.cpp libsais/liblibsais.a -lsdsl -ldivsufsort -ldivsufsort64 -DINFO_SPACE -o k2bp_info.x;

k2bp_comp_info:
	g++ -I ~/local-sdsl-lite/include -L ~/local-sdsl-lite/lib k2bp_compr_info.cpp libsais/liblibsais.a -lsdsl -ldivsufsort -ldivsufsort64 -DINFO_SPACE -o k2bp_compr_info.x;

k2bp_mult_1:
	g++ -I ~/local-sdsl-lite/include -L ~/local-sdsl-lite/lib k2bp_mult_1.cpp libsais/liblibsais.a -lsdsl -ldivsufsort -ldivsufsort64 -O3 -DNDEBUG -m64 -o k2bp_mult_1.x;

k2bp_mult_comp_1:
	g++ -I ~/local-sdsl-lite/include -L ~/local-sdsl-lite/lib k2bp_compr_mult.cpp libsais/liblibsais.a -lsdsl -ldivsufsort -ldivsufsort64 -O3 -DNDEBUG -m64 -o k2bp_compr_mult.x;

k2bp_rand_test:
	g++ -I ~/local-sdsl-lite/include -L ~/local-sdsl-lite/lib test.cpp libsais/liblibsais.a -lsdsl -ldivsufsort -ldivsufsort64 -m64 -O3 -DNDEBUG -o test.x;
