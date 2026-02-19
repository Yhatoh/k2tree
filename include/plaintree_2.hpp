#ifndef __PLAIN_TREE__
#define __PLAIN_TREE__

// std includes
#include <cstdlib>
#include <vector>

// sdsl includes
#include <sdsl/int_vector.hpp>

// local includes
#include "bit_vector.hpp"

using namespace sdsl;
using namespace std;

struct plain_tree {
  block_vector tree;
  block_vector l;

  uint8_t height_tree;
  uint64_t msize;
  uint64_t rmsize;
  uint64_t m;

  plain_tree() {}
  
  inline void destroy() {
    tree.destroy();
    l.destroy();
  }

  inline void binsum(plain_tree &B, plain_tree &C) {
    // trivial cases
    if(B.tree.size() == 2 && tree.size() == 2) {
      C.tree.push_back(1);
      C.tree.push_back(0);
      C.height_tree = height_tree;
      C.msize = msize;
      C.rmsize = rmsize;
      C.m = m;
      return;
    }

    if(B.tree.size() == 2) {
      C.tree.swap(tree);
      C.l.swap(l);
      C.height_tree = height_tree;
      C.msize = msize;
      C.rmsize = rmsize;
      C.m = m;
      return;
    }

    if(tree.size() == 2) {
      C.tree.swap(B.tree);
      C.l.swap(B.l);
      C.height_tree = B.height_tree;
      C.msize = B.msize;
      C.rmsize = B.rmsize;
      C.m = B.m;
      return;
    }

    uint8_t curr_depth = 0;
    tree.init_scan(); B.tree.init_scan();
    l.init_scan(); B.l.init_scan();

    C.height_tree = height_tree;
    C.msize = msize;
    C.rmsize = rmsize;
    C.m = m;
    return;
  }
};

#endif
