#include "k2tree_bp_sdsl_idems.hpp"
//#include "k2tree_bp.hpp"

template< typename T, typename T2 >
ostream& operator<<(ostream& os, const pair< T, T2 > &p) {
  os << "(" << p.first << "," << p.second << ")";
  return os;
}

template< typename T >
ostream& operator<<(ostream& os, const vector< T > &vec) {
  os << "[";
  for(uint64_t i = 0; i < vec.size(); i++) {
    os << vec[i];
    if(i != vec.size() - 1)
      os << ", ";
  }
  os << "]";
  return os;
}

int main() {

  vector< vector< uint64_t > > matrix = {
    {0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
  
  vector< pair< uint64_t, uint64_t > > ones;
  for(uint64_t i = 0; i < matrix.size(); i++) {
    for(uint64_t j = 0; j < matrix.size(); j++) {
      if(matrix[i][j]) ones.push_back({i, j});
    }
  }

  k2tree_bp_sdsl<> k2tree(ones);
  cout << k2tree << "\n";
  {
  auto check = k2tree.get_pos_ones();
  assert(check.size() == ones.size());
  sort(check.begin(), check.end());
  sort(ones.begin(), ones.end());
  for(uint64_t i = 0; i < ones.size(); i++) {
    assert(check[i] == ones[i]);
  }
  }

  k2tree_bp_sdsl_idems<2,
    sd_vector<>, rank_support_sd<1>,
    sd_vector<>, rank_support_sd<1>, rank_support_sd<0>,
    select_support_sd<1>, select_support_sd<0>> k2tree_idem(k2tree);
  cout << k2tree_idem << "\n";

  auto check = k2tree_idem.get_pos_ones();
  assert(check.size() == ones.size());
  sort(check.begin(), check.end());
  sort(ones.begin(), ones.end());
  cout << check << endl;
  cout << ones << endl;
  for(uint64_t i = 0; i < ones.size(); i++) {
    assert(check[i] == ones[i]);
  }
  return 0;
}
