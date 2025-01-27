// std includes
#include <iostream>

// local includes
#include "k2tree_bp_sdsl.hpp"
#include "k2tree_bp_sdsl_idems.hpp"

int main(int argc, char** argv) {
  if(argc <= 1) {
    std::cerr << "At least one arguments:" << endl;
    std::cerr << "  ./k2bp_build.x <path file matrix> <size of squared matrix> <amount of ones>" << endl;
    exit(1);
  }

  for(uint64_t i = 1; i < argc; i += 3) {
    cerr << "Getting parameters..." << endl;

    std::string k2_path = argv[i];

    ifstream k2_file;
    k2_file.open(k2_path);

    if(!k2_file.is_open()) {
      cerr << "Error opening file. Check if the file exists or the path is writed correctly" << endl;
      exit(1);
    }

    cerr << "Reading k2tree..." << endl;

    k2tree_bp_sdsl<2> k2tree;
    k2tree.load(k2_file);

    auto ones = k2tree.get_pos_ones();

    cerr << "Building k2tree..." << endl;


    cerr << "Checking if it is correct..." << endl;

    auto check = k2tree.get_pos_ones();
    sort(check.begin(), check.end());

    k2tree_bp_sdsl_idems<2,
      sd_vector<>, rank_support_sd<1>,
      sd_vector<>, rank_support_sd<1>, rank_support_sd<0>,
      select_support_sd<1>, select_support_sd<0>> k2tree_idem(k2tree);
    assert(check == ones);

    cerr << "Writing file..." << endl;

    ofstream k2_file_idem;
    k2_file.open(k2_path + ".k2bpi");
    k2tree_idem.write(k2_file);
    k2_file.close();

  }
  return 0;
}
