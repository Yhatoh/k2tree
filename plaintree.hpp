#ifndef __PLAIN_TREE__
#define __PLAIN_TREE__

// std includes
#include <vector>

// sdsl includes
#include <sdsl/int_vector.hpp>

using namespace sdsl;
using namespace std;

struct plain_tree {
  vector< uint8_t > tree;
  vector< uint8_t > l;

  uint64_t height_tree;
  uint64_t msize;
  uint64_t rmsize;
  uint64_t m;


  inline void add_one(vector< uint64_t > &bv, uint64_t &pos_to_add) {
    bv.push_back(pos_to_add++);
  }

  inline void add_zero(vector< uint64_t > &bv, uint64_t &pos_to_add) {
    pos_to_add++;
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
      swap(C.tree, tree);
      swap(C.l, l);
      C.height_tree = height_tree;
      C.msize = msize;
      C.rmsize = rmsize;
      C.m = m;
      return;
    }

    if(tree.size() == 2) {
      swap(C.tree, B.tree);
      swap(C.l, B.l);
      C.height_tree = B.height_tree;
      C.msize = B.msize;
      C.rmsize = B.rmsize;
      C.m = B.m;
      return;
    }
    uint64_t curr_bit_tree = 0;
    vector< uint64_t > bits_tree;

#ifdef DEBUG
    cout << "Starting algorithm" << endl;
#endif

    int64_t curr_depth = 0;
    while(A_tree < tree.size() && B_tree < B.tree.size()) {
      bool A_p = tree[A_tree];
      bool B_p = tree[B_tree];
      if(A_p && B_p && curr_depth < height_tree) {
#ifdef DEBUG
        cout << "Entering a subtree in both cases" << endl;
        cout << " curr depth: " << curr_depth << endl;
        cout << " pos tree A: " << A_tree << endl;
        cout << " pos L    A: " << A_L << endl;
        cout << " pos tree B: " << B_tree << endl;
        cout << " pos L    B: " << B_L << endl;
#endif
        A_tree++; B_tree++; curr_depth++;
        C.tree.push_back(1);
      } else if(A_p && B_p) {
#ifdef DEBUG
        cout << "Entering a subtree in both cases and last level" << endl;
        cout << " curr depth: " << curr_depth << endl;
        cout << " pos tree A: " << A_tree << endl;
        cout << " pos L    A: " << A_L << endl;
        cout << " pos tree B: " << B_tree << endl;
        cout << " pos L    B: " << B_L << endl;
#endif
        C.tree.push_back(1);
        C.l.push_back(l[A_L] | B.l[B_L]);
        A_L++;
        B_L++;
        A_tree++;
        B_tree++;
        curr_depth++;
      } else if(A_p && !B_p) {
#ifdef DEBUG
        cout << "Copying subtree A" << endl;
        cout << " curr depth: " << curr_depth << endl;
        cout << " pos tree A: " << A_tree << endl;
        cout << " pos L    A: " << A_L << endl;
        cout << " pos tree B: " << B_tree << endl;
        cout << " pos L    B: " << B_L << endl;
#endif
        // end tree A
        uint64_t counter = 1;
        while(counter > 0) {
          if(tree[A_tree] && curr_depth < height_tree) {
            C.tree.push_back(1);
            counter++;
            curr_depth++;
          } else if(tree[A_tree]) {
            C.tree.push_back(1);
            C.l.push_back(l[A_L]);
            A_L++;
            counter++;
            curr_depth++;
          } else {
            C.tree.push_back(0);
            counter--;
            curr_depth--;
          }
          A_tree++;
        }

        // end tree B
        B_tree++;
      } else if(!A_p && B_p) {
#ifdef DEBUG
        cout << "Copying subtree B" << endl;
        cout << " curr depth: " << curr_depth << endl;
        cout << " pos tree A: " << A_tree << endl;
        cout << " pos L    A: " << A_L << endl;
        cout << " pos tree B: " << B_tree << endl;
        cout << " pos L    B: " << B_L << endl;
#endif
        // end tree A
        A_tree++;

        // end tree B
        uint64_t counter = 1;
        while(counter > 0) {
          if(B.tree[B_tree] && curr_depth < height_tree) {
            C.tree.push_back(1);
            counter++;
            curr_depth++;
          } else if(B.tree[B_tree]) {
            C.l.push_back(B.l[B_L]);
            C.tree.push_back(1);
            B_L++;
            counter++;
            curr_depth++;
          } else {
            C.tree.push_back(0);
            counter--;
            curr_depth--;
          }
          B_tree++;
        }
      } else {
#ifdef DEBUG
        cout << "Both submatrices 0" << endl;
        cout << " curr depth: " << curr_depth << endl;
        cout << " pos tree A: " << A_tree << endl;
        cout << " pos L    A: " << A_L << endl;
        cout << " pos tree B: " << B_tree << endl;
        cout << " pos L    B: " << B_L << endl;
#endif
        C.tree.push_back(0);
        A_tree++; B_tree++;
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
