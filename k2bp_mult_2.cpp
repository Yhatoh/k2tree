// std includes
#include <iostream>

// local includes
#include "k2tree_bp_sdsl.hpp"

int main(int argc, char** argv) {
  if(argc <= 1) {
    std::cerr << "At least two arguments" << endl;
    std::cerr << "  ./k2bp_build.x <path file matrix1.k2bp> ... <path file matrixn.k2bp" << endl;
    exit(1);
  }

  std::string k2_path = argv[1];

  std::ifstream k2_file;
  k2_file.open(k2_path);

  if(!k2_file.is_open()) {
    cerr << "Error opening file. Check if the file exists or the path is writed correctly" << endl;
    exit(1);
  }

  k2tree_bp_sdsl<2> result;
  result.load(k2_file);
  k2_file.close();

  for(uint64_t i = 2; i < argc; i++) {
    std::string k2_path = argv[i];

    std::ifstream k2_file;
    k2_file.open(k2_path);

    if(!k2_file.is_open()) {
      cerr << "Error opening file. Check if the file exists or the path is writed correctly" << endl;
      exit(1);
    }

    k2tree_bp_sdsl<2> k2tree;
    k2tree.load(k2_file);
    k2_file.close();

    result = result * k2tree;
  }

  std::stringstream name_file;
  name_file << k2_path << ".mul2";

  std::ofstream result_file;
  result_file.open(name_file.str());
  result.write(result_file);
  result_file.close();

  return 0;
}
