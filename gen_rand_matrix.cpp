#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

int main(int argv, char* argc[]) {
  if(argv < 3 || argv > 5) {
    std::cerr << "Expect at least 3 arguments at maximum 5 argument:" << std::endl;
    std::cerr << "  ./gen_rand_matrix.x <size> <density> <amount> (path) (seed)" << std::endl;
    exit(1);
  }

  uint64_t seed = 42;
  if(argv == 3) {
    seed = atoi(argc[2]);
  }
  std::srand(seed);

  uint64_t size = atoi(argc[1]);
  double p = std::stod(argc[2]);
  uint64_t amount = atoi(argc[3]);
  uint64_t m = size * p;
  std::string path = argc[4];
  for(uint64_t i = 0; i < amount; i++) {
    std::vector< uint8_t > matrix(size * size, 0);
    std::fill(matrix.begin(), matrix.begin() + m, 1);
    std::random_shuffle(matrix.begin(), matrix.end());

    std::stringstream path_save;
    path_save << path;
    path_save << "/rand_matrix." << argc[2] << "." << i << ".txt";
    std::ofstream save_matrix(path_save.str());

    for(uint64_t i = 0; i < size; i++) {
      for(uint64_t j = 0; j < size; j++) {
        if(matrix[i * size + j])
          save_matrix << i << " " << j << "\n";
      }
    }
    save_matrix.close();
  }
}
