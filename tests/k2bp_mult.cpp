// std includes
#include <iostream>

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

  while((c=getopt(argc, argv, "b")) != -1) {
    switch(c) {
      case 'b':
        type = 1; break;
      case '?':
        cerr << "Unkown option: " << optarg << std::endl;
        exit(1);
    }
  }

  optind -= 1;
  if(argc - optind != 3) usage_and_exit(argv[0]);
  argv += optind; argc -= optind;

  std::string k2_1_path = argv[1];
  std::string k2_2_path = argv[2];

  std::ifstream k2_1_file;
  k2_1_file.open(k2_1_path);

  if(!k2_1_file.is_open()) {
    cerr << "Error opening file. Check if the file exists or the path is writed correctly" << endl;
    exit(1);
  }

  std::ifstream k2_2_file;
  k2_2_file.open(k2_2_path);

  if(!k2_2_file.is_open()) {
    cerr << "Error opening file. Check if the file exists or the path is writed correctly" << endl;
    exit(2);
  }

  if(type == 0) {
    k2_bp<2, bit_vector> m1;
    m1.load(k2_1_file);
    k2_1_file.close();

    k2_bp<2, bit_vector> m2;
    m2.load(k2_2_file);
    k2_2_file.close();

    k2_bp<2, bit_vector> m3;
    m1.mul(m2, m3);
    

    std::stringstream name_file;
    name_file << k2_1_path << ".mul";

    std::ofstream result_file;
    result_file.open(name_file.str());
    m3.write(result_file);
    result_file.close();
  } else {
//    k2_bp<2, rrr_vector<127>> m1;
//    m1.load(k2_1_file);
//    k2_1_file.close();
//
//    k2_bp<2, rrr_vector<127>> m2;
//    m2.load(k2_2_file);
//    k2_2_file.close();
//
//    k2_bp<2, rrr_vector<127>> m3;
//    m1.mul(m2, m3);
//    
//
//    std::stringstream name_file;
//    name_file << k2_1_path << ".mul";
//
//    std::ofstream result_file;
//    result_file.open(name_file.str());
//    m3.write(result_file);
//    result_file.close();
  }
  return 0;
}

void usage_and_exit(char* name) {
    cerr << "Usage:\n\t  " << name <<  " [options] infile1 infile2 \n\n";
    cerr << "Options:\n";
    cerr << "\t-b        use rrr_vector<127> to store leaves (def. bit_vector) [bit_vector more space, faster; rrr_vector less space, slower]\n";    
    cerr << "Multiply two compressed matrices stored in infile1 and infile2\n\n";
    exit(1);
}
