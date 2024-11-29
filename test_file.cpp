#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <string>

#include "k2tree_bp_sdsl.hpp"
#include "k2tree_bp_sdsl_idems.hpp"

using namespace std;

int main(int argc, char** argv) {

  string n = argv[1];
  for(uint64_t i = 2; i < argc; i += 1) {
    ifstream file;
    file.open(argv[i]);

    uint64_t amount = 0;
    vector< pair< uint64_t, uint64_t > > ones;
    
    cout << "Reading file " << argv[i] << endl;
    assert(file.is_open());
    while(!file.eof()) {
      uint64_t a, b;

      file >> a >> b;

      if(file.eof()) break;
      ones.push_back({a, b});
      amount++;
    }

    sort(ones.begin(), ones.end());
    file.close();

    cout << amount << endl;

    k2tree_bp_sdsl<2> k2tree(ones, stoi(n));

    cout << k2tree.size_in_bits() << " " << (double) k2tree.size_in_bits() / amount << " " << (double) k2tree.size_in_bits() / k2tree.nodes() << endl;
//    auto ret = k2tree.get_pos_ones();
//
//    sort(ret.begin(), ret.end());
//    assert(ret.size() == ones.size());

    k2tree_bp_sdsl_idems<2, sd_vector<>, rank_support_sd<>, rank_support_sd<0>, select_support_sd<>, select_support_sd<0>> k2tree_idems(k2tree);
    cout << k2tree_idems.size_in_bits() << " " << (double) k2tree_idems.size_in_bits() / amount << " " << (double) k2tree_idems.size_in_bits() / k2tree.nodes() << endl;

//    auto ret = k2tree.get_pos_ones();
//    sort(ret.begin(), ret.end());
//    assert(ret.size() == ones.size());


  }
  return 0;
}
