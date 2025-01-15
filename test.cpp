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
Randomer genmatrix(256, 256, 49);
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
  assert(check.size() == ones.size());
  sort(ones.begin(), ones.end());
  sort(check.begin(), check.end());
  for(uint64_t i = 0; i < ones.size(); i++) {
    assert(check[i] == ones[i]);
  }
  //cout << k2tree << endl;

  k2tree_bp_sdsl_idems<2,
    sd_vector<>, rank_support_sd<1>,
    sd_vector<>, rank_support_sd<1>, rank_support_sd<0>,
    select_support_sd<1>, select_support_sd<0>> k2tree_idem(k2tree);
  //cout << k2tree_idem << endl;
  auto check2 = k2tree_idem.get_pos_ones();

  assert(check2.size() == ones.size());

  sort(check2.begin(), check2.end());
  for(uint64_t i = 0; i < ones.size(); i++) {
    assert(check2[i] == ones[i]);
  }

  return true;
}

bool test_gen_matrices(uint64_t n, uint64_t m) {
  cout << "Generating a matrix of size " << n << "x"  << m << endl;
  vector< pair< uint64_t, uint64_t> > ones = gen_ones_matrix(n, m);

  cout << "Generating k2 tree" << endl;
  k2tree_bp_sdsl<2> k2tree(ones);
  k2tree_bp_sdsl_idems<2,
    sd_vector<>, rank_support_sd<1>,
    sd_vector<>, rank_support_sd<1>, rank_support_sd<0>,
    select_support_sd<1>, select_support_sd<0>> k2tree_idem(k2tree);

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

bool test_union_algorithm(uint64_t n, uint64_t m) {
  cout << "Generating a matrix of size " << n << "x"  << m << endl;
  vector< pair< uint64_t, uint64_t> > ones = gen_ones_matrix(n, m);
  vector< pair< uint64_t, uint64_t> > ones2 = gen_ones_matrix(n, m);

  cout << "Generating A" << endl;
  k2tree_bp_sdsl<2> A(ones);
  cout << "Generating B" << endl;
  k2tree_bp_sdsl<2> B(ones2);

  cout << "Getting ones from k2 tree" << endl;
  auto ones_A = A.get_pos_ones();
  auto ones_B = B.get_pos_ones();
  set<pair<uint64_t, uint64_t>> union_;
  for(auto p : ones_A) union_.insert(p);
  for(auto p : ones_B) union_.insert(p);

  vector< pair< uint64_t, uint64_t > > ones_union;
  for(auto p : union_) ones_union.push_back(p);

  auto C = A | B;

  auto check = C.get_pos_ones();

  assert(check.size() == ones_union.size());

  sort(check.begin(), check.end());
  sort(ones_union.begin(), ones_union.end());
  for(uint64_t i = 0; i < ones_union.size(); i++) {
    assert(check[i] == ones_union[i]);
  }

  return true;
}

bool test_multi_algorithm(uint64_t n, uint64_t m) {
  cout << "Generating a matrix of size " << n << "x"  << m << endl;
  vector< pair< uint64_t, uint64_t> > ones = gen_ones_matrix(n, m);
  vector< pair< uint64_t, uint64_t> > ones2 = gen_ones_matrix(n, m);
  vector< pair< uint64_t, uint64_t> > expected;

  cout << "Brute force multiplication" << endl;
  {
    vector< vector< uint64_t > > mA(n, vector<uint64_t>(m, 0));
    for(auto p : ones) mA[p.first][p.second] = 1;

    vector< vector< uint64_t > > mB(n, vector<uint64_t>(m, 0));
    for(auto p : ones2) mB[p.first][p.second] = 1;

    vector< vector< uint64_t > > mC(n, vector< uint64_t >(m, 0));

    for(uint64_t i = 0; i < n; i++) {
      for(uint64_t j = 0; j < m; j++) {
        for(uint64_t k = 0; k < n; k++) {
          mC[i][j] |= mA[i][k] & mB[k][j];
        }
      }
    }


    for(uint64_t i = 0; i < n; i++) {
      for(uint64_t j = 0; j < m; j++) {
        if(mC[i][j]) expected.push_back({i, j});
      }
    }
  }

  cout << "Generating A" << endl;
  k2tree_bp_sdsl<2> A(ones);
  cout << "Generating B" << endl;
  k2tree_bp_sdsl<2> B(ones2);

  cout << "C = A * B" << endl;
  auto C = A * B;

  auto check = C.get_pos_ones();

  assert(check.size() == expected.size());

  sort(check.begin(), check.end());
  sort(expected.begin(), expected.end());
  for(uint64_t i = 0; i < expected.size(); i++) {
    assert(check[i] == expected[i]);
  }

  return true;
}

bool test_union_algorithm_tree_comp(uint64_t n, uint64_t m) {
  cout << "Generating a matrix of size " << n << "x"  << m << endl;
  vector< pair< uint64_t, uint64_t> > ones = gen_ones_matrix(n, m);
  vector< pair< uint64_t, uint64_t> > ones2 = gen_ones_matrix(n, m);

  cout << "Generating A" << endl;
  k2tree_bp_sdsl<2> A(ones);
  k2tree_bp_sdsl_idems<2,
    sd_vector<>, rank_support_sd<1>,
    sd_vector<>, rank_support_sd<1>, rank_support_sd<0>,
    select_support_sd<1>, select_support_sd<0>> A_idem(A);
  cout << "Generating B" << endl;
  k2tree_bp_sdsl<2> B(ones2);
  k2tree_bp_sdsl_idems<2,
    sd_vector<>, rank_support_sd<1>,
    sd_vector<>, rank_support_sd<1>, rank_support_sd<0>,
    select_support_sd<1>, select_support_sd<0>> B_idem(B);

  cout << "Getting ones from A" << endl;
  auto ones_A = A_idem.get_pos_ones();
  cout << "Getting ones from B" << endl;
  auto ones_B = B_idem.get_pos_ones();
  cout << "Getting ones from A | B (bruteforce)" << endl;
  set<pair<uint64_t, uint64_t>> union_;
  for(auto p : ones_A) union_.insert(p);
  for(auto p : ones_B) union_.insert(p);

  vector< pair< uint64_t, uint64_t > > ones_union;
  for(auto p : union_) ones_union.push_back(p);

  cout << "A | B" << endl;
  auto C = A_idem | B_idem;
  cout << "Compressing C" << endl;
  k2tree_bp_sdsl_idems<2,
    sd_vector<>, rank_support_sd<1>,
    sd_vector<>, rank_support_sd<1>, rank_support_sd<0>,
    select_support_sd<1>, select_support_sd<0>> C_idem(C);

  auto check = C_idem.get_pos_ones();

  assert(check.size() == ones_union.size());

  sort(check.begin(), check.end());
  sort(ones_union.begin(), ones_union.end());
  for(uint64_t i = 0; i < ones_union.size(); i++) {
    assert(check[i] == ones_union[i]);
  }

  return true;
}


int main() {
//  cout << "Testing pow 2 square matrices" << endl;
//  for(uint64_t t = 0; t < 50; t++) {
//    cout << "Test " << t + 1 << endl;
//    test_pow_2_matrices(pow2matrix());
//    cout << " Passed!" << endl;
//  }
//  cout << "Testing general matrices" << endl;
//  for(uint64_t t = 0; t < 50; t++) {
//    cout << "Test " << t + 1 << endl;
//    test_gen_matrices(genmatrix(), genmatrix());
//    cout << "Passed!" << endl;
//  }
//
//  cout << "Testing Union Algorithm" << endl;
//  for(uint64_t t = 0; t < 50; t++) {
//    cout << "Test " << t + 1 << endl;
//    test_union_algorithm(genmatrix(), genmatrix());
//    cout << "Passed!" << endl;
//  }

  cout << "Testing Union Algorithm tree compression" << endl;
  for(uint64_t t = 0; t < 50; t++) {
    cout << "Test " << t + 1 << endl;
    test_union_algorithm_tree_comp(genmatrix(), genmatrix());
    cout << "Passed!" << endl;
  }
//
//  cout << "Testing Multiply Algorithm" << endl;
//  for(uint64_t t = 0; t < 50; t++) {
//    cout << "Test " << t + 1 << endl;
//    test_multi_algorithm(genmatrix(), genmatrix());
//    cout << "Passed!" << endl;
//  }
  return 0;
}
