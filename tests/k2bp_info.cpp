// std includes
#include <iostream>

// local includes
#include "k2_bp.hpp"

void usage_and_exit(char* name);

int main(int argc, char** argv) {
  extern char *optarg;
  extern int optind, opterr, optopt;

  int type = 0;
  int c;

  while((c=getopt(argc, argv, "bt:cs:")) != -1) {
    switch(c) {
      case 'b':
        type = 1; break;
      case '?':
        cerr << "Unkown option: " << optarg << std::endl;
        exit(1);
    }
  }

  optind -= 1;
  if(argc - optind != 2) usage_and_exit(argv[0]);
  argv += optind; argc -= optind;
  std::string matrix = argv[1];
  std::ifstream k2_file(matrix);
  
  if(type == 0) {
    k2_bp<2, bit_vector> k2tree;
    k2tree.load(k2_file);

    k2_file.close();


    cout << "Information k2tree " << matrix << endl;
    cout << " Size Matrix: " << k2tree.size_matrix() << " Amount of 1's: " << k2tree.size() << endl;
    cout << " Amount of nodes: " << k2tree.nodes() << endl;
    uint64_t bits = k2tree.size_in_bits();
    cout << " Bits    : " << bits << endl;
    cout << " Bits/1's: " << (double) bits / k2tree.size() << endl;
    cout << " Bits/n  : " << (double) bits / k2tree.nodes() << endl;
  } else { 
    k2_bp<2, rrr_vector<127>> k2tree;
    k2tree.load(k2_file);

    k2_file.close();


    cout << "Information k2tree " << matrix << endl;
    cout << " Size Matrix: " << k2tree.size_matrix() << " Amount of 1's: " << k2tree.size() << endl;
    cout << " Amount of nodes: " << k2tree.nodes() << endl;
    uint64_t bits = k2tree.size_in_bits();
    cout << " Bits    : " << bits << endl;
    cout << " Bits/1's: " << (double) bits / k2tree.size() << endl;
    cout << " Bits/n  : " << (double) bits / k2tree.nodes() << endl;
  }

  return 0;
}

void usage_and_exit(char* name) {
    cerr << "Usage:\n\t  " << name <<  " [options] filename \n\n";
    cerr << "Options:\n";
    cerr << "\t-b        use rrr_vector<127> to store leaves (def. bit_vector) [bit_vector more space, faster; rrr_vector less space, slower]\n";    
    cerr << "Show info about compress matrix in filename\n\n";
    exit(1);
}
