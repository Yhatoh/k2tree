// std includes
#include <iostream>
#include <utility>
#include <set>

// local includes
#include "k2tree_bp_sdsl.hpp"
//#include "k2tree_bp_sdsl_intL.hpp"

int main(int argc, char** argv) {
  if(argc < 3) {
    std::cerr << "At least three arguments:" << endl;
    std::cerr << "  ./k2bp_build.x <path file matrix> <size of squared matrix> <amount of ones>" << endl;
    exit(1);
  }

  if((argc - 1) % 3) {
    std::cerr << "The arguments should be of group of 3:" << endl;
    std::cerr << "  ./k2bp_build.x <matrix1> <size1> <#ones1> ... <matrixn> <sizen> <#onesn>" << endl;
    exit(1);
  }

  for(uint64_t i = 1; i < argc; i += 3) {
    cerr << "Getting parameters..." << endl;

    std::string matrix = argv[i];
    uint64_t size = std::stoi(argv[i + 1]);
    uint64_t m = std::stoi(argv[i + 2]);

    vector< pair< uint64_t, uint64_t > > ones;

    ifstream ones_txt;
    ones_txt.open(matrix);

    if(!ones_txt.is_open()) {
      cerr << "Error opening file. Check if the file exists or the path is writed correctly" << endl;
      exit(1);
    }

    cerr << "Reading ones..." << endl;
    for(uint64_t i = 0; i < m; i++) {
      uint64_t x, y;
      ones_txt >> x >> y;
      ones.push_back({x, y});
    }
    sort(ones.begin(), ones.end());

    ones_txt.close();
    cout << ones.size() << endl;

    cerr << "Building k2tree..." << endl;

    k2tree_bp_sdsl<2, rrr_vector<127>> k2tree(ones, size);

    cerr << "Checking if it is correct..." << endl;

    auto check = k2tree.get_pos_ones();
    sort(check.begin(), check.end());

    assert(check == ones);

    cerr << "Writing file..." << endl;

    ofstream k2_file;
    k2_file.open(matrix + ".k2bp");
    k2tree.write(k2_file);
    k2_file.close();

  }
  return 0;
}
