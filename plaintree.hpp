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


  void add_one(vector< uint64_t > &bv, uint64_t &pos_to_add) {
    bv.push_back(pos_to_add++);
  }

  void add_zero(vector< uint64_t > &bv, uint64_t &pos_to_add) {
    pos_to_add++;
  }

  void binsum(plain_tree &B, plain_tree &C) {
    uint64_t pa, pb;
    uint64_t pLa, pLb;

    pa = pb = pLa = pLb = 0;

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
    while(pa < tree.size() && pb < B.tree.size()) {
      if(tree[pa] && B.tree[pb] && curr_depth < height_tree) {
#ifdef DEBUG
        cout << "Entering a subtree in both cases" << endl;
        cout << " curr depth: " << curr_depth << endl;
        cout << " pos tree A: " << pa << endl;
        cout << " pos L    A: " << pLa << endl;
        cout << " pos tree B: " << pb << endl;
        cout << " pos L    B: " << pLb << endl;
#endif
        pa++; pb++; curr_depth++;
        add_one(bits_tree, curr_bit_tree);
      } else if(tree[pa] && B.tree[pb]) {
#ifdef DEBUG
        cout << "Entering a subtree in both cases and last level" << endl;
        cout << " curr depth: " << curr_depth << endl;
        cout << " pos tree A: " << pa << endl;
        cout << " pos L    A: " << pLa << endl;
        cout << " pos tree B: " << pb << endl;
        cout << " pos L    B: " << pLb << endl;
#endif
        add_one(bits_tree, curr_bit_tree);
        for(uint64_t i = 0; i < 4; i++) {
          if((pLa < l.size() && l[pLa]) || 
              (pLb < B.l.size() && B.l[pLb])) {
            add_one(bits_L, curr_bit_L);
          } else {
            add_zero(bits_L, curr_bit_L);
          }
          pLa++; pLb++;
        }
        pa++;
        pb++;
        curr_depth++;
      } else if(tree[pa] && !B.tree[pb]) {
#ifdef DEBUG
        cout << "Copying subtree A" << endl;
        cout << " curr depth: " << curr_depth << endl;
        cout << " pos tree A: " << pa << endl;
        cout << " pos L    A: " << pLa << endl;
        cout << " pos tree B: " << pb << endl;
        cout << " pos L    B: " << pLb << endl;
#endif
        // end tree A
        uint64_t counter = 1;
        while(counter > 0) {
          if(tree[pa] && curr_depth < height_tree) {
            add_one(bits_tree, curr_bit_tree);
            counter++;
            curr_depth++;
          } else if(tree[pa]) {
            add_one(bits_tree, curr_bit_tree);
            for(uint64_t i = 0; i < 4; i++) {
              if(l[pLa]) add_one(bits_L, curr_bit_L);
              else add_zero(bits_L, curr_bit_L);
              pLa++;
            }
            counter++;
            curr_depth++;
          } else {
            add_zero(bits_tree, curr_bit_tree);
            counter--;
            curr_depth--;
          }
          pa++;
        }

        // end tree B
        pb++;
      } else if(!tree[pa] && B.tree[pb]) {
#ifdef DEBUG
        cout << "Copying subtree B" << endl;
        cout << " curr depth: " << curr_depth << endl;
        cout << " pos tree A: " << pa << endl;
        cout << " pos L    A: " << pLa << endl;
        cout << " pos tree B: " << pb << endl;
        cout << " pos L    B: " << pLb << endl;
#endif
        // end tree A
        pa++;

        // end tree B
        uint64_t counter = 1;
        while(counter > 0) {
          if(B.tree[pb] && curr_depth < height_tree) {
            add_one(bits_tree, curr_bit_tree);
            counter++;
            curr_depth++;
          } else if(B.tree[pb]) {
            add_one(bits_tree, curr_bit_tree);
            for(uint64_t i = 0; i < 4; i++) {
              if(B.l[pLb]) add_one(bits_L, curr_bit_L);
              else add_zero(bits_L, curr_bit_L);
              pLb++;
            }
            counter++;
            curr_depth++;
          } else {
            add_zero(bits_tree, curr_bit_tree);
            counter--;
            curr_depth--;
          }
          pb++;
        }
      } else {
#ifdef DEBUG
        cout << "Both submatrices 0" << endl;
        cout << " curr depth: " << curr_depth << endl;
        cout << " pos tree A: " << pa << endl;
        cout << " pos L    A: " << pLa << endl;
        cout << " pos tree B: " << pb << endl;
        cout << " pos L    B: " << pLb << endl;
#endif
        add_zero(bits_tree, curr_bit_tree);
        pa++; pb++;
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
