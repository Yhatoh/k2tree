// c++ includes
#include <algorithm>
#include <iostream>

// sdsl includes
#include <sdsl/int_vector.hpp>

// local includes
#include "bit_vector.hpp"

sdsl::bit_vector random_bit_vector(size_t n, size_t seed = 42) {
  sdsl::bit_vector bits(n, 0);
  sdsl::util::set_random_bits(bits, seed);
  return bits;
}

#define N 100000
#define NN 100
int main() {
  {
    std::cout << "testing access" << std::endl;
    auto bv_test = random_bit_vector(N);

    bvector bv_min;
    for(size_t i = 0; i < N; i++) {
      bv_min.push_back(bv_test[i]);
    }

    for(size_t i = 0; i < N; i++) {
      assert(bv_min[i] == bv_test[i]);
    }
    std::cout << "test passed" << std::endl;
  }
  {
    std::cout << "testing read int" << std::endl;
    auto bv_test = random_bit_vector(N);

    bvector bv_min;
    for(size_t i = 0; i < N; i++) {
      bv_min.push_back(bv_test[i]);
    }

    for(size_t i = 0; i < N - 64; i++) {
      uint64_t len = (std::rand() % 64) + 1;
      auto res = bv_min.read_int(i, len);
      assert(res == bv_test.get_int(i, len));
    }
    std::cout << "test passed" << std::endl;
  }
  {
    std::cout << "testing concat" << std::endl;

    bvector bv_min;
    std::vector< bool > test_bool;
    for(size_t j = 0; j < NN; j++) {
      auto bv_test = random_bit_vector(N);
      bvector bv_min2;
      for(size_t i = 0; i < N; i++) {
        bv_min2.push_back(bv_test[i]);
      }

      size_t x = (std::rand() % bv_test.size());
      size_t len = (std::rand() % (std::min((uint64_t) 64, (uint64_t) bv_test.size() - x)) + 1);
      for(size_t i = 0; i < len; i++) {
        test_bool.push_back(bv_test[x + i]);
      }
      bv_min.concat(bv_min2, x, x + len);
    }

    for(size_t i = 0; i < test_bool.size(); i++) {
      assert(bv_min[i] == test_bool[i]);
    }
    std::cout << "test passed" << std::endl;
  }
}
