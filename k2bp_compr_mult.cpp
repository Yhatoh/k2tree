// std includes
#include <iostream>

// local includes
#include "k2tree_bp_sdsl_idems.hpp"

int main(int argc, char** argv) {
  if(argc <= 1 || argc > 3) {
    std::cerr << "Two arguments expected" << endl;
    std::cerr << "  ./k2bp_build.x <path file matrix1.k2bp> <path file matrix2.k2bp>" << endl;
    exit(1);
  }

  std::string k2_1_path = argv[1];
  std::string k2_2_path = argv[2];

  std::ifstream k2_1_file;
  k2_1_file.open(k2_1_path);

  if(!k2_1_file.is_open()) {
    cerr << "Error opening file. Check if the file exists or the path is writed correctly" << endl;
    exit(1);
  }

  k2tree_bp_sdsl_idems<2, rrr_vector<127>,
    sd_vector<>, rank_support_sd<1>,
    sd_vector<>, rank_support_sd<1>, rank_support_sd<0>,
    select_support_sd<1>, select_support_sd<0>> m1;
  m1.load(k2_1_file);
  k2_1_file.close();

  std::ifstream k2_2_file;
  k2_2_file.open(k2_2_path);

  if(!k2_2_file.is_open()) {
    cerr << "Error opening file. Check if the file exists or the path is writed correctly" << endl;
    exit(2);
  }

  k2tree_bp_sdsl_idems<2, rrr_vector<127>,
    sd_vector<>, rank_support_sd<1>,
    sd_vector<>, rank_support_sd<1>, rank_support_sd<0>,
    select_support_sd<1>, select_support_sd<0>> m2;
  m2.load(k2_2_file);
  k2_2_file.close();

  plain_tree result;
  m1.mul(m2, result);
  k2tree_bp_sdsl<2, rrr_vector<127>> m3(result);
  

  std::stringstream name_file;
  name_file << k2_1_path << ".mulrrri";

  std::ofstream result_file;
  result_file.open(name_file.str());
  m3.write(result_file);
  result_file.close();

  return 0;
}
