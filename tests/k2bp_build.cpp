// std includes
#include <getopt.h>
#include <iostream>
#include <utility>
#include <set>

// local includes
#include "k2_bp.hpp"

void usage_and_exit(char* argv);

int main(int argc, char** argv) {
  extern char *optarg;
  extern int optind, opterr, optopt;

  int type = 0;
  float p = -1LL;
  int64_t size = -1LL;
  bool check = 0;
  int c;

  while((c=getopt(argc, argv, "bt:cs:")) != -1) {
    switch(c) {
      case 'b':
        type = 1; break;
      case 't':
        p = std::stof(optarg); break;
      case 'c':
        check = true; break;
      case 's':
        size = std::stoll(optarg); break;
      case '?':
        cerr << "Unkown option: " << optarg << std::endl;
        exit(1);
    }
  }

  optind -= 1;
  if(argc - optind != 2) usage_and_exit(argv[0]);
  argv += optind; argc -= optind;
  // checking values of args
  if(p != -1LL) {
    if(p < 0) {
      cerr << "percentage of threshold (-t " << p << ") should not be negative" << std::endl;
      exit(1);
    }
  }
  if(size != -1) {
    if(size <= 0) {
      cerr << "size of matrix (-s " << size << ") should at least 1 (obviously)" << std::endl;
      exit(1);
    }
  }


  std::string matrix = argv[1];

  vector< pair< uint64_t, uint64_t > > ones;

  ifstream ones_txt;
  ones_txt.open(matrix);

  if(!ones_txt.is_open()) {
    cerr << "Error opening file. Check if the file exists or the path is writed correctly" << endl;
    exit(1);
  }

  uint64_t x, y;
  uint64_t max_x = 0, max_y = 0;
  while(ones_txt >> x >> y) {
    ones.push_back({x, y});
    if(size == -1LL) {
      max_x = std::max(x, max_x);
      max_y = std::max(y, max_y);
    }
  }
  if(size == -1LL) {
    size = std::max(max_x, max_y) + 1;
  }
  //sort(ones.begin(), ones.end());

  ones_txt.close();

  if(type == 0) {
    k2_bp<2, bit_vector> k2tree(ones, size);
    if(p != -1) {
      k2tree.add_child_info(std::sqrt(k2tree.nodes()) * p);
    }

    std::vector< std::pair< uint64_t, uint64_t > > check;
    k2tree.get_pos_ones(check);
    sort(check.begin(), check.end());

    assert(check == ones);

    ofstream k2_file;
    k2_file.open(matrix + ".k2bp");
    k2tree.write(k2_file);
    k2_file.close();
  } else {
    k2_bp<2, rrr_vector<127>> k2tree(ones, size);
    if(p != -1) {
      k2tree.add_child_info(std::sqrt(k2tree.nodes()) * p);
    }

    std::vector< std::pair< uint64_t, uint64_t > > check;
    k2tree.get_pos_ones(check);
    sort(check.begin(), check.end());

    assert(check == ones);

    ofstream k2_file;
    k2_file.open(matrix + ".k2bp");
    k2tree.write(k2_file);
    k2_file.close();
  }
  return 0;
}

void usage_and_exit(char* name) {
    cerr << "Usage:\n\t  " << name << " [options] filename \n\n";
    cerr << "Options:\n";
    cerr << "\t-b        use rrr_vector<127> to store leaves (def. bit_vector) [bit_vector more space, faster; rrr_vector less space, slower]\n";    
    cerr << "\t-c        compressed->decompress->check\n";
    cerr << "\t-s S      matrix actual size (def. largest index + 1)\n";
    cerr << "\t-t p      use p * sqrt(S) as threshold for subtree information (def. don't add subtree information)\n";    
    cerr << "Compress filename.k2bp\n\n";
    exit(1);
}
