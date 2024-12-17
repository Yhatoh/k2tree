compile_test:
	g++ -I ~/local-sdsl-lite/include -L ~/local-sdsl-lite/lib test_file.cpp libsais/liblibsais.a -lsdsl -ldivsufsort -ldivsufsort64 -o out

compile_test_local:
	g++ -I ~/include -L ~/lib test_file.cpp libsais/liblibsais.a -lsdsl -ldivsufsort -ldivsufsort64 -o out
