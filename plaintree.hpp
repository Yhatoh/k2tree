#ifndef __PLAIN_TREE__
#define __PLAIN_TREE__

// std includes
#include <vector>

// sdsl includes
#include <sdsl/int_vector.hpp>

using namespace sdsl;
using namespace std;

struct plain_tree {
  sdsl::bit_vector tree;
  sdsl::bit_vector l;

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
      C.tree = bit_vector(2, 0);
      C.tree[0] = 1;
      C.height_tree = height_tree;
      C.msize = msize;
      C.rmsize = rmsize;
      C.m = m;
      return;
    }

    if(B.tree.size() == 2) {
      C.tree = tree;
      C.l = l;
      C.height_tree = height_tree;
      C.msize = msize;
      C.rmsize = rmsize;
      C.m = m;
      return;
    }

    if(tree.size() == 2) {
      C.tree = B.tree;
      C.l = B.l;
      C.height_tree = B.height_tree;
      C.msize = B.msize;
      C.rmsize = B.rmsize;
      C.m = B.m;
      return;
    }
    uint64_t curr_bit_tree = 0;
    vector< uint64_t > bits_tree;

    uint64_t curr_bit_L = 0;
    vector< uint64_t > bits_L;

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
        add_one(bits_tree, curr_bit_tree);
      } else if(A_p && B_p) {
#ifdef DEBUG
        cout << "Entering a subtree in both cases and last level" << endl;
        cout << " curr depth: " << curr_depth << endl;
        cout << " pos tree A: " << A_tree << endl;
        cout << " pos L    A: " << A_L << endl;
        cout << " pos tree B: " << B_tree << endl;
        cout << " pos L    B: " << B_L << endl;
#endif
        add_one(bits_tree, curr_bit_tree);
        for(uint64_t i = 0; i < 4; i++) {
          if((A_L < l.size() && l[A_L]) || 
              (B_L < B.l.size() && B.l[B_L])) {
            add_one(bits_L, curr_bit_L);
          } else {
            add_zero(bits_L, curr_bit_L);
          }
          A_L++; B_L++;
        }
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
            add_one(bits_tree, curr_bit_tree);
            counter++;
            curr_depth++;
          } else if(tree[A_tree]) {
            add_one(bits_tree, curr_bit_tree);
            for(uint64_t i = 0; i < 4; i++) {
              if(l[A_L]) add_one(bits_L, curr_bit_L);
              else add_zero(bits_L, curr_bit_L);
              A_L++;
            }
            counter++;
            curr_depth++;
          } else {
            add_zero(bits_tree, curr_bit_tree);
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
            add_one(bits_tree, curr_bit_tree);
            counter++;
            curr_depth++;
          } else if(B.tree[B_tree]) {
            add_one(bits_tree, curr_bit_tree);
            if(B.l[B_L]) add_one(bits_L, curr_bit_L);
            else add_zero(bits_L, curr_bit_L);

            if(B.l[B_L + 1]) add_one(bits_L, curr_bit_L);
            else add_zero(bits_L, curr_bit_L);

            if(B.l[B_L + 2]) add_one(bits_L, curr_bit_L);
            else add_zero(bits_L, curr_bit_L);

            if(B.l[B_L + 3]) add_one(bits_L, curr_bit_L);
            else add_zero(bits_L, curr_bit_L);

            B_L += 4;
            counter++;
            curr_depth++;
          } else {
            add_zero(bits_tree, curr_bit_tree);
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
        add_zero(bits_tree, curr_bit_tree);
        A_tree++; B_tree++;
        curr_depth--;
      }
    }


    C.tree = bit_vector(curr_bit_tree, 0);
    for(auto bit : bits_tree) C.tree[bit] = 1;

    C.l = bit_vector(curr_bit_L, 0);
    for(auto bit : bits_L) C.l[bit] = 1;

    C.height_tree = height_tree;
    C.msize = msize;
    C.rmsize = rmsize;
    C.m = m;
    return;
  }
};

#endif
