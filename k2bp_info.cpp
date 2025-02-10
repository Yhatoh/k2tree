// std includes
#include <iostream>

// local includes
#include "k2tree_bp_sdsl.hpp"
//#include "k2tree_bp_sdsl_intL.hpp"

int main(int argc, char** argv) {
  if(argc <= 1) {
    std::cerr << "At least one argument:" << endl;
    std::cerr << "  ./k2bp_info.x <path file k2tree bp>" << endl;
    exit(1);
  }

  for(uint64_t i = 1; i < argc; i += 3) {
    std::string k2_path = argv[i];

    ifstream k2_file;
    k2_file.open(k2_path);

    if(!k2_file.is_open()) {
      cerr << "Error opening file. Check if the file exists or the path is writed correctly" << endl;
      exit(1);
    }

    k2tree_bp_sdsl<2, sd_vector<>> k2tree;
    k2tree.load(k2_file);

    k2_file.close();


    cout << "Information k2tree " << k2_path << endl;
    cout << " Size Matrix: " << k2tree.size_matrix() << " Amount of 1's: " << k2tree.size() << endl;
    cout << " Amount of nodes: " << k2tree.nodes() << endl;
    uint64_t bits = k2tree.size_in_bits();
    cout << " Bits    : " << bits << endl;
    cout << " Bits/1's: " << (double) bits / k2tree.size() << endl;
    cout << " Bits/n  : " << (double) bits / k2tree.nodes() << endl;

  }
  return 0;
}
