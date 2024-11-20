#ifndef K2_TREE_BP
#define K2_TREE_BP

// std includes
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <stack>
#include <utility>
#include <vector>
#include <tuple>

// local includes
#include "util.hpp"

using namespace std;

// k2-tree
// parameters:
//   * k * k: amount of children per node
template< uint64_t k = 2 >
class k2tree_bp {
  private:
    uint64_t height_tree;

    vector< uint8_t > tree; // k2tree
    uint8_t last_bit_t;

    vector< uint8_t > l; // real values
    uint8_t last_bit_l;

    uint64_t msize;
    uint64_t rmsize;
    uint64_t m;

    bool has_at_least_one(vector< vector< uint64_t > > &adj_list, vector< vector< uint64_t > ::iterator > &pointers, uint64_t init_x, uint64_t init_y, uint64_t subm_size) {
      for(uint64_t x = init_x; x < init_x + subm_size; x++) {
        if(pointers[x] != adj_list[x].end() && *(pointers[x]) < init_y + subm_size) {
          return true;
        }
      }
      return false;
    }

    void add_one(vector< uint8_t > &bv, uint64_t &pos_to_add) {
      bv.back() |= 1 << pos_to_add;
      pos_to_add++;
      if(pos_to_add == 8) {
        bv.push_back(0);
        pos_to_add = 0;
      }
    }

    void add_zero(vector< uint8_t > &bv, uint64_t &pos_to_add) {
      pos_to_add++;
      if(pos_to_add == 8) {
        bv.push_back(0);
        pos_to_add = 0;
      }
    }

  public:
    k2tree_bp(vector< pair< uint64_t, uint64_t > > &ones, uint64_t n = -1) { 
      m = ones.size();

      if(n == -1) {
        // minimum size of a matriz, max index + 1
        n = 0;
        for(const auto& one : ones) n = max(one.first, max(one.second, n)); 
      }

      rmsize = n;
      height_tree = ceil_log2(n);
      msize = (1 << height_tree);

#ifdef DEBUG
      cerr << "Real Size Matrix: " << rmsize << "x" << rmsize << endl;
      cerr << "Size Matrix: " << msize << "x" << msize << endl;
      cerr << "Height Tree: " << height_tree << endl;
#endif // DEBUG
      
      sort(ones.begin(), ones.end());

#ifdef DEBUG
      cerr << "Generating adjacency list..." << endl;
#endif // DEBUG
      vector< vector< uint64_t > > adj_list(msize);
      for(const auto& one : ones) adj_list[one.first].push_back(one.second);

#ifdef DEBUG
      cerr << "Generating pointers..." << endl;
#endif // DEBUG
      vector< vector< uint64_t >::iterator > pointers(msize);
      for(uint64_t i = 0; i < msize; i++) pointers[i] = adj_list[i].begin();

      // first i will do it asuming k = 2
      // then i will generalize
#ifdef DEBUG
      cerr << "Initialize recursion..." << endl;
#endif // DEBUG
      stack< tuple< uint64_t, uint64_t, uint64_t, bool > > recursion;
      recursion.push(make_tuple(msize / 2, msize / 2, msize / 2, false));
      recursion.push(make_tuple(msize / 2, msize / 2, 0, false));
      recursion.push(make_tuple(msize / 2, 0, msize / 2, false));
      recursion.push(make_tuple(msize / 2, 0, 0, false));

      uint64_t pos_to_add = 1;
      uint64_t pos_to_add_l = 0;

#ifdef DEBUG
      cerr << "Adding ( of root..." << endl;
#endif // DEBUG
      tree.push_back(1);
      l.push_back(0);

      while(!recursion.empty()) {
        auto [subm_size, init_x, init_y, flag] = recursion.top();
#ifdef DEBUG
        cerr << "Current Sub Matrix Size: " << subm_size << " init x: " << init_x << " init y: " << init_y << " visited: " << flag << endl;
#endif // DEBUG
        recursion.pop();

        if(subm_size == k && has_at_least_one(adj_list, pointers, init_x, init_y, subm_size)) { // submatrix of size k^2
          
          for(uint64_t x = init_x; x < init_x + subm_size; x++) {
            for(uint64_t y = init_y; y < init_y + subm_size; y++) {
              if(pos_to_add_l == 8) {
                l.push_back(0);
                pos_to_add_l = 0;
              }
              if(pointers[x] != adj_list[x].end() && *(pointers[x]) == y) {
                l.back() |= 1 << pos_to_add_l;
                pointers[x]++;
                pos_to_add_l++;
              } else {
                pos_to_add_l++;
              }
            }
          }

#ifdef DEBUG
          cerr << "Adding (..." << endl;
#endif // DEBUG
          add_one(tree, pos_to_add);
#ifdef DEBUG
          cerr << "Adding (..." << endl;
#endif // DEBUG
          add_one(tree, pos_to_add);
#ifdef DEBUG
          cerr << "Adding )..." << endl;
#endif // DEBUG
          add_zero(tree, pos_to_add);
#ifdef DEBUG
          cerr << "Adding )..." << endl;
#endif // DEBUG
          add_zero(tree, pos_to_add);
        } else if(subm_size == k) {
#ifdef DEBUG
          cerr << "Adding (..." << endl;
#endif // DEBUG
          add_one(tree, pos_to_add);
#ifdef DEBUG
          cerr << "Adding )..." << endl;
#endif // DEBUG
          add_zero(tree, pos_to_add);
        } else {

          if(flag) {
#ifdef DEBUG
          cerr << "Adding )..." << endl;
#endif // DEBUG
            add_zero(tree, pos_to_add);
            continue;
          } else {
            recursion.push(make_tuple(subm_size, init_x, init_y, true));
          }

#ifdef DEBUG
          cerr << "Adding (..." << endl;
#endif // DEBUG
          add_one(tree, pos_to_add);
          
          if(has_at_least_one(adj_list, pointers, init_x, init_y, subm_size)) {
            recursion.push(make_tuple(subm_size / 2, init_x + subm_size / 2, init_y + subm_size / 2, false));
            recursion.push(make_tuple(subm_size / 2, init_x + subm_size / 2, init_y, false));
            recursion.push(make_tuple(subm_size / 2, init_x, init_y + subm_size / 2, false));
            recursion.push(make_tuple(subm_size / 2, init_x, init_y, false));
          } else {
#ifdef DEBUG
          cerr << "Adding )..." << endl;
#endif // DEBUG
            add_zero(tree, pos_to_add);
            recursion.pop();
          }
        }
      }
#ifdef DEBUG
      cerr << "Adding ) of root ..." << endl;
#endif // DEBUG
      pos_to_add++;
      if(pos_to_add == 8) {
        tree.push_back(0);
        pos_to_add = 0;
      }
      recursion.pop();

      last_bit_t = pos_to_add;
      last_bit_l = pos_to_add_l;
    }

    vector< pair< uint64_t, uint64_t > > get_pos_ones() {
      stack< tuple< uint8_t, uint64_t, uint64_t > > child_visit;

      uint64_t r, c;
      r = c = 0;
      child_visit.push({0, r, c});

      uint64_t to_read_l = 0;
      vector< pair< uint64_t, uint64_t > > ret;

      for(uint64_t i = 0; i < tree.size(); i++) {
        for(uint8_t j = (i == 0 ? 1 : 0); j < (tree.size() - 1 == i ? last_bit_t : 8); j++) {
#ifdef DEBUG
          cerr << "Total bits " << (tree.size() - 1) * 8 + last_bit_t + 1 << endl;
          cerr << "Reading " << i * 8 + j << " bit" << endl;
#endif // DEBUG
          if(tree[i] & (1 << j)) {
            auto [vis, r_, c_] = child_visit.top();
            r = r_ + (vis / k) * (1 << (height_tree - child_visit.size()));
            c = c_ + (vis % k) * (1 << (height_tree - child_visit.size()));
            child_visit.push({0, r, c});
#ifdef DEBUG
            cerr << "Level " << child_visit.size() << endl;
            cerr << "Current row: " << r << " col: " << c << endl;
#endif // DEBUG
          } else {
#ifdef DEBUG
            cerr << "End of subtree" << endl;
#endif // DEBUG
            if(child_visit.size() == height_tree + 1) {
#ifdef DEBUG
              cerr << "Last level, read real values" << endl;
              cerr << "L size " << (l.size() - 1) * 8 + last_bit_l << " ";
              cerr << "Current bit " << to_read_l << endl;
#endif // DEBUG
              for(uint8_t i = 0; i < k * k; i++) {
                if(l[to_read_l / 8] & (1 << (to_read_l % 8))) {
                  auto [vis, r_, c_] = child_visit.top();
                  ret.push_back({r_ + i / k, c_ + i % k});
                }
                to_read_l++;
              }
#ifdef DEBUG
              cerr << "Finishing reading real values" << endl;
#endif // DEBUG
            }
            child_visit.pop();
            // means we finish to read the complete tree
            if(child_visit.size() != 0) {
              auto [vis, r_, c_] = child_visit.top();
              child_visit.pop();
              child_visit.push({vis + 1, r_, c_});
            }
          }
        }
      }
      return ret;
    }

    friend ostream& operator<<(ostream& os, const k2tree_bp<k> &k2tree) {
      cout << "Height Tree: " << k2tree.height_tree << endl;
      cout << "Tree: ";
      for(uint64_t i = 0; i < k2tree.tree.size(); i++) {
        for(uint8_t j = 0; j < (k2tree.tree.size() - 1 == i ? k2tree.last_bit_t : 8); j++) {
          cout << ((k2tree.tree[i] & (1 << j)) ? "(" : ")");
        }
      }
      cout << endl;
      cout << "L: ";
      for(uint64_t i = 0; i < k2tree.l.size(); i++) {
        for(uint8_t j = 0; j < (k2tree.l.size() - 1 == i ? k2tree.last_bit_l : 8); j++) {
          if(j % 4 == 0 && !(i == 0 && j == 0)) cout << " ";
          cout << ((k2tree.l[i] & (1 << j)) ? "1" : "0");
        }
      }
      return os;
    }
};
#endif // !K2_TREE_BP
