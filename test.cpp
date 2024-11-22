#include "k2tree_bp_sdsl.hpp"
#include "k2tree_bp_sdsl_idems.hpp"
//#include "k2tree_bp.hpp"


#include <algorithm>
#include <random>

class Randomer {
  // random seed by default
  std::mt19937 gen_;
  std::uniform_int_distribution< size_t > dist_;

  public:
    
    Randomer(size_t min, size_t max, unsigned int seed = std::random_device{}()) 
      : gen_{seed}, dist_(min, max) {}

    // if you want predictable numbers 
    void SetSeed(unsigned int seed) {
      gen_.seed(seed);
    }

    size_t operator()() {
      return dist_(gen_);
    }
    
};

Randomer pow2matrix(2, 10, 49);
Randomer genmatrix(1, 10000, 49);
Randomer zerone(0, 1, 49);

vector< pair< uint64_t, uint64_t > > gen_ones_matrix(uint64_t n, uint64_t m) {
    vector< pair< uint64_t, uint64_t > > ones;
    for (uint64_t i = 0; i < n; ++i) {
        for (uint64_t j = 0; j < m; ++j) {
            if(zerone()) {
              ones.push_back({i, j});
            }
        }
    }
    return ones;
}

bool test_pow_2_matrices(uint64_t p) {
  cout << "Generating a matrix of size " << (1 << p) << "x"  << (1 << p) << endl;
  vector< pair< uint64_t, uint64_t> > ones = gen_ones_matrix(1 << p, 1 << p);

  cout << "Generating k2 tree" << endl;
  k2tree_bp_sdsl<2> k2tree(ones);

  cout << "Getting ones from k2 tree" << endl;

  auto check = k2tree.get_pos_ones();

  k2tree_bp_sdsl_idems<2, sd_vector<>, rank_support_sd<>, rank_support_sd<0>, select_support_sd<>, select_support_sd<0>> k2tree_idem(k2tree);

  auto check2 = k2tree_idem.get_pos_ones();

  assert(check.size() == ones.size());
  assert(check2.size() == ones.size());

  sort(check.begin(), check.end());
  sort(check2.begin(), check2.end());
  sort(ones.begin(), ones.end());
  for(uint64_t i = 0; i < ones.size(); i++) {
    assert(check[i] == ones[i]);
    assert(check2[i] == ones[i]);
  }

  return true;
}

bool test_gen_matrices(uint64_t n, uint64_t m) {
  cout << "Generating a matrix of size " << n << "x"  << m << endl;
  vector< pair< uint64_t, uint64_t> > ones = gen_ones_matrix(n, m);

  cout << "Generating k2 tree" << endl;
  k2tree_bp_sdsl<2> k2tree(ones);
  k2tree_bp_sdsl_idems<2, sd_vector<>, rank_support_sd<>, rank_support_sd<0>, select_support_sd<>, select_support_sd<0>> k2tree_idem(k2tree);

  cout << "Getting ones from k2 tree" << endl;
  auto check = k2tree.get_pos_ones();
  auto check2 = k2tree_idem.get_pos_ones();

  assert(check.size() == ones.size());
  assert(check2.size() == ones.size());

  sort(check.begin(), check.end());
  sort(check2.begin(), check2.end());
  sort(ones.begin(), ones.end());
  for(uint64_t i = 0; i < ones.size(); i++) {
    assert(check[i] == ones[i]);
    assert(check2[i] == ones[i]);
  }

  return true;
}


int main() {
  cout << "Testing pow 2 square matrices" << endl;
  for(uint64_t t = 0; t < 50; t++) {
    cout << "Test " << t + 1 << endl;
    test_pow_2_matrices(pow2matrix());
    cout << " Passed!" << endl;
  }
  cout << "Testing general matrices" << endl;
  for(uint64_t t = 0; t < 50; t++) {
    cout << "Test " << t + 1 << endl;
    test_gen_matrices(genmatrix(), genmatrix());
    cout << "Passed!" << endl;
  }
  return 0;
}
