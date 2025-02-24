#include "k2tree_bp_sdsl.hpp"
#include "k2tree_bp_sdsl_idems.hpp"

#include <algorithm>
#include <fstream>
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
Randomer genmatrix(8024, 8024, 49);
Randomer zerone(0, 10000, 49);
vector< pair< uint64_t, uint64_t > > gen_ones_matrix(uint64_t n, uint64_t m) {
    vector< pair< uint64_t, uint64_t > > ones;
    for (uint64_t i = 0; i < n; ++i) {
        for (uint64_t j = 0; j < m; ++j) {
            if(zerone() > 9990) {
              ones.push_back({i, j});
            }
        }
    }
    cout << ones.size() << endl;
    return ones;
}

bool test_pow_2_matrices(uint64_t p) {
  cout << "Generating a matrix of size " << (1 << p) << "x"  << (1 << p) << endl;
  vector< pair< uint64_t, uint64_t> > ones = gen_ones_matrix(1 << p, 1 << p);

  cout << "Generating k2 tree" << endl;
  k2tree_bp_sdsl<2, bit_vector> k2tree(ones);

  cout << "Getting ones from k2 tree" << endl;

  auto check = k2tree.get_pos_ones();
  assert(check.size() == ones.size());
  sort(ones.begin(), ones.end());
  sort(check.begin(), check.end());
  for(uint64_t i = 0; i < ones.size(); i++) {
    assert(check[i] == ones[i]);
  }

  k2tree_bp_sdsl_idems<2, bit_vector,
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
  k2tree_bp_sdsl<2, bit_vector> k2tree(ones);
  k2tree_bp_sdsl_idems<2, bit_vector,
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

bool test_multi_algorithm(uint64_t n, uint64_t m) {
  cout << "Generating a matrix of size " << n << "x"  << m << endl;
  vector< pair< uint64_t, uint64_t> > ones = gen_ones_matrix(n, m);
  vector< pair< uint64_t, uint64_t> > ones2 = gen_ones_matrix(n, m);
#ifdef DEBUG
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
#endif

  cout << "Generating A" << endl;
  k2tree_bp_sdsl<2, rrr_vector<127>> A(ones);
  cout << "Generating B" << endl;
  k2tree_bp_sdsl<2, rrr_vector<127>> B(ones2);

  cout << "C = A * B" << endl;
  plain_tree aux_C;
  A.mul(B, aux_C);
  k2tree_bp_sdsl<2, rrr_vector<127>> C(aux_C);

#ifdef DEBUG
  auto check = C.get_pos_ones();

  assert(check.size() == expected.size());

  sort(check.begin(), check.end());
  sort(expected.begin(), expected.end());
  for(uint64_t i = 0; i < expected.size(); i++) {
    assert(check[i] == expected[i]);
  }
#endif

  return true;
}

bool test_multi_algorithm_tree_comp(uint64_t n, uint64_t m) {
  cout << "Generating a matrix of size " << n << "x"  << m << endl;
  vector< pair< uint64_t, uint64_t> > ones = gen_ones_matrix(n, m);
  vector< pair< uint64_t, uint64_t> > ones2 = gen_ones_matrix(n, m);
#ifdef DEBUG
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
#endif

  cout << "Generating A" << endl;
  k2tree_bp_sdsl<2, bit_vector> A(ones);
  k2tree_bp_sdsl_idems<2, bit_vector,
    sd_vector<>, rank_support_sd<1>,
    sd_vector<>, rank_support_sd<1>, rank_support_sd<0>,
    select_support_sd<1>, select_support_sd<0>> A_idem(A);

  cout << "Generating B" << endl;
  k2tree_bp_sdsl<2, bit_vector> B(ones2);
  k2tree_bp_sdsl_idems<2, bit_vector,
    sd_vector<>, rank_support_sd<1>,
    sd_vector<>, rank_support_sd<1>, rank_support_sd<0>,
    select_support_sd<1>, select_support_sd<0>> B_idem(A);

  cout << "C = A * B" << endl;
  plain_tree aux_C;
  A.mul(B, aux_C);
  k2tree_bp_sdsl<2, bit_vector> C(aux_C);
  cout << "Compressing C" << endl;
  k2tree_bp_sdsl_idems<2, bit_vector,
    sd_vector<>, rank_support_sd<1>,
    sd_vector<>, rank_support_sd<1>, rank_support_sd<0>,
    select_support_sd<1>, select_support_sd<0>> C_idem(C);

#ifdef DEBUG

  auto check = C_idem.get_pos_ones();

  assert(check.size() == expected.size());

  sort(check.begin(), check.end());
  sort(expected.begin(), expected.end());
  for(uint64_t i = 0; i < expected.size(); i++) {
    assert(check[i] == expected[i]);
  }
#endif // DEBUG

  return true;
}

bool test_write_load(uint64_t n, uint64_t m) {
  cout << "Generating a matrix of size " << n << "x"  << m << endl;
  vector< pair< uint64_t, uint64_t> > ones = gen_ones_matrix(n, m);

  cout << "Generating k2 tree" << endl;
  k2tree_bp_sdsl<2, bit_vector> k2tree(ones);
  {
    cout << "Writing on file" << endl;
    ofstream file_k2_write;
    file_k2_write.open("matrix.k2");
    k2tree.write(file_k2_write);
    file_k2_write.close();
    cout << "Reading from file" << endl;

    ifstream file_k2_read;
    file_k2_read.open("matrix.k2");
    k2tree_bp_sdsl<2, bit_vector> k2tree_load;
    k2tree_load.load(file_k2_read);
    file_k2_read.close();

    cout << "Check ones k2" << endl;

    auto check = k2tree_load.get_pos_ones();
    sort(check.begin(), check.end());
    sort(ones.begin(), ones.end());
    for(uint64_t i = 0; i < ones.size(); i++) {
      assert(check[i] == ones[i]);
    }
  }
  cout << "Compressing k2 tree" << endl;
  k2tree_bp_sdsl_idems<2, bit_vector,
    sd_vector<>, rank_support_sd<1>,
    sd_vector<>, rank_support_sd<1>, rank_support_sd<0>,
    select_support_sd<1>, select_support_sd<0>> k2tree_idem(k2tree);
  {
    cout << "Writing on file" << endl;
    ofstream file_k2_write;
    file_k2_write.open("matrix.k2_idem");
    k2tree_idem.write(file_k2_write);
    file_k2_write.close();
    cout << "Reading from file" << endl;

    ifstream file_k2_read;
    file_k2_read.open("matrix.k2");
    k2tree_bp_sdsl<2, bit_vector> k2tree_load;
    k2tree_load.load(file_k2_read);
    file_k2_read.close();

    cout << "Check ones k2idem" << endl;

    auto check = k2tree_load.get_pos_ones();
    sort(check.begin(), check.end());
    sort(ones.begin(), ones.end());
    for(uint64_t i = 0; i < ones.size(); i++) {
      assert(check[i] == ones[i]);
    }
  }

  return true;
}

int main(int argc, char *argv[]) {
  bool f_compr = 0;
  bool f_multi = 0;
  bool f_multi_compr = 0;
  bool f_write_and_load = 0;
  if(argc == 1) {
    f_compr = f_multi = f_multi_compr = f_write_and_load = 1;
  } else {
    int c;
    while ((c=getopt(argc, argv, "tmcw")) != -1) {
      switch (c) {
        case 't':
          f_compr = true; break;
        case 'm':
          f_multi = true; break;
        case 'c':
          f_multi_compr = true; break;
        case 'w':
          f_write_and_load = true; break;
        case '?':
          fprintf(stderr,"Unknown option: %c\n", optopt);
          exit(1);
      }
    }
  }
  if(f_compr) {
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
  }

  if(f_multi) {
    cout << "Testing Multiply Algorithm" << endl;
    for(uint64_t t = 0; t < 1; t++) {
      cout << "Test " << t + 1 << endl;
      test_multi_algorithm(genmatrix(), genmatrix());
      cout << "Passed!" << endl;
    }
  }

  if(f_multi_compr) {
    cout << "Testing Multiply Algorithm tree compression" << endl;
    for(uint64_t t = 0; t < 1; t++) {
      cout << "Test " << t + 1 << endl;
      test_multi_algorithm_tree_comp(genmatrix(), genmatrix());
      cout << "Passed!" << endl;
    }
  }
  
  if(f_write_and_load) {
    cout << "Testing writing and load" << endl;
    for(uint64_t t = 0; t < 50; t++) {
      cout << "Test " << t + 1 << endl;
      test_write_load(genmatrix(), genmatrix());
      cout << "Passed!" << endl;
    }
  }
  return 0;
}
