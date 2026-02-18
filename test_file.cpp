#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <string>

#include "k2tree_bp_sdsl.hpp"
#include "k2tree_bp_sdsl_idems.hpp"

#include <sdsl/rrr_vector.hpp>

using namespace std;

int main(int argc, char** argv) {

  for(uint64_t i = 1; i < argc; i += 2) {
    string n = argv[i];
    ifstream file;
    file.open(argv[i + 1]);

    uint64_t amount = 0;
    vector< pair< uint64_t, uint64_t > > ones;
    
    cout << "Reading file " << argv[i + 1] << endl;
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

    k2tree_bp_sdsl<2, sd_vector<>> k2tree(ones, stoi(n));

    uint64_t amount_bits = k2tree.size_in_bits();
    cout << amount_bits << " " << (double) amount_bits / amount << " " << (double) amount_bits / k2tree.nodes() << endl;
    auto ret = k2tree.get_pos_ones();
//
    //sort(ret.begin(), ret.end());
    //assert(ret.size() == ones.size());

//    k2tree_bp_sdsl_idems<2,
//                         rrr_vector<127>, rank_support_rrr<1, 127>,
//                         rrr_vector<127>, rank_support_rrr<1, 127>, rank_support_rrr<0, 127>,
//                                          select_support_rrr<1, 127>, select_support_rrr<0, 127>> k2tree_idems(k2tree);
    k2tree_bp_sdsl_idems<2, sd_vector<>,
                         sd_vector<>, rank_support_sd<>,
                         sd_vector<>, rank_support_sd<>, rank_support_sd<0>,
                                          select_support_sd<1>, select_support_sd<0>> k2tree_idems(k2tree);

    uint64_t amount_bits_ = k2tree_idems.size_in_bits();
    cout << amount_bits_ << " " << (double) amount_bits_ / amount << " " << (double) amount_bits_ / k2tree.nodes() << endl;

    //ret = k2tree.get_pos_ones();
    //sort(ret.begin(), ret.end());
    //assert(ret.size() == ones.size());


  }
  return 0;
}
