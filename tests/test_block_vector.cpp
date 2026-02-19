// c++ includes
#include <iostream>
#include <random>
#include <vector>

// local includes
#include "bit_vector.hpp"

std::vector<bool> random_bit_vector(size_t n, size_t seed = 42) {
  static std::mt19937 gen(seed);        // Mersenne Twister RNG
  static std::bernoulli_distribution d(0.5); // 50/50 chance for 0 or 1

  std::vector<bool> bits(n);
  for (size_t i = 0; i < n; ++i) {
    bits[i] = d(gen);
  }
  return bits;
}

#define N 100000
#define NN 100
int main() {

  {
    std::cout << "testing push_back" << std::endl;
    auto bv_test = random_bit_vector(N);

    block_vector bv;

    for(auto x : bv_test) {
      bv.push_back(x);
    }

    bv.init_scan();

    for(size_t i = 0; i < bv_test.size(); i++) {
      auto obtain = bv.read();
      bv.move();
      assert(obtain == bv_test[i]);
    }

    std::cout << "all test passed" << std::endl;
  }
  {
    std::cout << "testing concat" << std::endl;

    std::vector< bool > bv_big_test;
    std::vector< bool > bv_test = random_bit_vector(N);

    block_vector bv_big;
    for(auto x : bv_test) bv_big.push_back(x);
    bv_big_test.insert(bv_big_test.end(), bv_test.begin(), bv_test.end());

    for(size_t i = 0; i < NN; i++) {
      block_vector bv_small;
      std::vector< bool > bv_test = random_bit_vector(N);

      for(auto x : bv_test) bv_small.push_back(x);
      bv_big_test.insert(bv_big_test.end(), bv_test.begin(), bv_test.end());

      bv_big.concat(bv_small);
    }

    bv_big.init_scan();
    for(size_t i = 0; i < bv_big_test.size(); i++) {
      assert(bv_big.read() == bv_big_test[i]);
      bv_big.move();
    }

    std::cout << "all test passed" << std::endl;
  }
  {
    std::cout << "testing copy" << std::endl;

    std::vector< bool > bv_big_test;
    std::vector< bool > bv_test = random_bit_vector(N);

    block_vector bv_big;
    for(auto x : bv_test) bv_big.push_back(x);
    bv_big_test.insert(bv_big_test.end(), bv_test.begin(), bv_test.end());

    for(size_t i = 0; i < NN; i++) {
      block_vector bv_small;
      std::vector< bool > bv_test = random_bit_vector(N);

      for(auto x : bv_test) bv_small.push_back(x);
      bv_big_test.insert(bv_big_test.end(), bv_test.begin(), bv_test.end());

      bv_big.concat(bv_small);
    }

    for(size_t i = 0; i < N; i++) {
      size_t x = 0 + (std::rand() % bv_big_test.size());
      size_t y = x + (std::rand() % (bv_big_test.size() - x));
      bv_big.init_scan();
      for(size_t i = 0; i < x; i++) bv_big.move();
      bv_big.save();

      for(size_t i = x; i < y; i++) bv_big.move();

      block_vector bv_range;
      bv_big.append_in(bv_range);

      bv_range.init_scan();
      for(size_t i = x; i <= y; i++) {
        assert(bv_range.read() == bv_big_test[i]);
        bv_range.move();
      }
    }

    std::cout << "all test passed" << std::endl;
  }
}
