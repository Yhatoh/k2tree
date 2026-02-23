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
  bvector tree;
  //bvector l;
  vector< uint8_t > l;

  uint8_t height_tree;
  uint64_t msize;
  uint64_t rmsize;
  uint64_t m;

  plain_tree() {}
  plain_tree(uint64_t n, uint64_t m) {
    tree.reserve(n);
    l.reserve(m);
  }

  void clear() {
    tree.clear();
    l.clear();
  }

  void reserve(uint64_t n, uint64_t m) {
    tree.reserve(n);
    l.reserve(m);
  }
  
  inline void destroy() {
    tree.destroy();
    //l.destroy();
    vector< uint8_t >().swap(l);
  }

  inline void binsum(plain_tree &B, plain_tree &C) {
    uint64_t A_tree, B_tree;
    uint64_t A_L, B_L;

    A_tree = B_tree = A_L = B_L = 0;

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
      swap(C.tree.bv, tree.bv);
      swap(C.tree.n, tree.n);
      swap(C.l, l);
      //swap(C.l.bv, l.bv);
      //swap(C.l.n, l.n);
      C.height_tree = height_tree;
      C.msize = msize;
      C.rmsize = rmsize;
      C.m = m;
      return;
    }

    if(tree.size() == 2) {
      swap(C.tree.bv, B.tree.bv);
      swap(C.tree.n, B.tree.n);
      swap(C.l, B.l);
      //swap(C.l.bv, B.l.bv);
      //swap(C.l.n, B.l.n);
      C.height_tree = B.height_tree;
      C.msize = B.msize;
      C.rmsize = B.rmsize;
      C.m = B.m;
      return;
    }

#ifdef DEBUG
    //cout << "Starting algorithm" << endl;
#endif

    uint8_t curr_depth = 0;
    for(; A_tree < tree.size() && B_tree < B.tree.size(); A_tree++, B_tree++) {
      uint64_t A_pos = A_tree;
      uint64_t B_pos = B_tree;
      uint64_t leaves = 0;
      for(; A_pos < tree.size() && B_pos < B.tree.size(); A_pos++, B_pos++) {
        if(tree[A_pos] != B.tree[B_pos]) break;
        curr_depth += (tree[A_pos] ? 1 : -1);
        //leaves += (curr_depth == height_tree + 1);
        if(curr_depth == height_tree + 1) {
          l[A_L + leaves] |= B.l[B_L + leaves];
          leaves++;
        }
      }

      // copy identical part and leaves
      C.tree.concat(tree, A_tree, A_pos);
      //for(uint32_t B_x = 0; B_x < leaves; B_x++) {
      //  l[A_L + B_x] |= B.l[B_L + B_x];
      //}
      C.l.insert(C.l.end(), l.begin() + A_L, l.begin() + A_L + leaves);

      A_tree = A_pos;
      B_tree = B_pos;
      A_L += leaves;
      B_L += leaves;

      // copy B
      if(A_tree < tree.size() && !tree[A_tree]) {
        B_pos = B_tree;
        uint8_t counter = curr_depth + 1;
        leaves = 0;
        B_pos++;
        leaves += (counter == height_tree + 1);
        for(; B_pos < B.tree.size() && counter >= curr_depth; B_pos++) {
          counter += (B.tree[B_pos] ? 1 : -1);
          leaves += (counter == height_tree + 1);
        }

        C.tree.concat(B.tree, B_tree, B_pos);
        //C.l.concat(B.l, B_L * 4, (B_L + leaves) * 4);
        C.l.insert(C.l.end(), B.l.begin() + B_L, B.l.begin() + B_L + leaves);

        B_tree = B_pos - 1;
        B_L += leaves;
        curr_depth--;
      }
      // copy A
      else if(B_tree < B.tree.size() && !B.tree[B_tree]) {
        A_pos = A_tree;
        uint8_t counter = curr_depth + 1;
        leaves = 0;
        A_pos++;
        leaves += (counter == height_tree + 1);
        for(; A_pos < tree.size() && counter >= curr_depth; A_pos++) {
          counter += (tree[A_pos] ? 1 : -1);
          leaves += (counter == height_tree + 1);
        }

        C.tree.concat(tree, A_tree, A_pos);
        C.l.insert(C.l.end(), l.begin() + A_L, l.begin() + A_L + leaves);

        A_tree = A_pos - 1;
        A_L += leaves;
        curr_depth--;
      }
    }

    C.height_tree = height_tree;
    C.msize = msize;
    C.rmsize = rmsize;
    C.m = m;
    return;
  }
};

#endif
