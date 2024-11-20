#ifndef K2_TREE_BP_SDSL
#define K2_TREE_BP_SDSL

// std includes
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <sdsl/construct.hpp>
#include <sdsl/io.hpp>
#include <sdsl/lcp_support_sada.hpp>
#include <stack>
#include <utility>
#include <vector>
#include <tuple>

// local includes
#include "util.hpp"

// sdsl includes
#include <sdsl/int_vector.hpp>
#include <sdsl/bp_support_sada.hpp>
#include <sdsl/util.hpp>
#include <sdsl/lcp.hpp>

using namespace std;
using namespace sdsl;

// k2-tree
// parameters:
//   * k * k: amount of children per node
template< uint64_t k = 2 >
class k2tree_bp_sdsl {
  public:
    uint64_t height_tree;

    bp_support_sada<> tree_support;
    bit_vector tree; // k2tree
    uint64_t last_bit_t; // universe

    bit_vector l; // real values
    uint64_t last_bit_l; // universe

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

    void add_one(vector< uint64_t > &bv, uint64_t &pos_to_add) {
      bv.push_back(pos_to_add++);
    }

    void add_zero(vector< uint64_t > &bv, uint64_t &pos_to_add) {
      pos_to_add++;
    }

  public:
    k2tree_bp_sdsl(vector< pair< uint64_t, uint64_t > > &ones, uint64_t n = -1) { 
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
      cout << "Real Size Matrix: " << rmsize << "x" << rmsize << endl;
      cout << "Size Matrix: " << msize << "x" << msize << endl;
      cout << "Height Tree: " << height_tree << endl;
#endif // DEBUG
      
      sort(ones.begin(), ones.end());

#ifdef DEBUG
      cout << "Generating adjacency list..." << endl;
#endif // DEBUG
      vector< vector< uint64_t > > adj_list(msize);
      for(const auto& one : ones) adj_list[one.first].push_back(one.second);

#ifdef DEBUG
      cout << "Generating pointers..." << endl;
#endif // DEBUG
      vector< vector< uint64_t >::iterator > pointers(msize);
      for(uint64_t i = 0; i < msize; i++) pointers[i] = adj_list[i].begin();

      // first i will do it asuming k = 2
      // then i will generalize
#ifdef DEBUG
      cout << "Initialize recursion..." << endl;
#endif // DEBUG
      stack< tuple< uint64_t, uint64_t, uint64_t, bool > > recursion;
      recursion.push(make_tuple(msize / 2, msize / 2, msize / 2, false));
      recursion.push(make_tuple(msize / 2, msize / 2, 0, false));
      recursion.push(make_tuple(msize / 2, 0, msize / 2, false));
      recursion.push(make_tuple(msize / 2, 0, 0, false));


#ifdef DEBUG
      cout << "Adding ( of root..." << endl;
#endif // DEBUG
      vector< uint64_t > bv_tree;
      vector< uint64_t > bv_l;

      bv_tree.push_back(0);

      uint64_t pos_to_add = 1;
      uint64_t pos_to_add_l = 0;

      while(!recursion.empty()) {
        auto [subm_size, init_x, init_y, flag] = recursion.top();
#ifdef DEBUG
        cout << "Current Sub Matrix Size: " << subm_size << " init x: " << init_x << " init y: " << init_y << " visited: " << flag << endl;
#endif // DEBUG
        recursion.pop();

        if(subm_size == k && has_at_least_one(adj_list, pointers, init_x, init_y, subm_size)) { // submatrix of size k^2
          
          for(uint64_t x = init_x; x < init_x + subm_size; x++) {
            for(uint64_t y = init_y; y < init_y + subm_size; y++) {
              if(pointers[x] != adj_list[x].end() && *(pointers[x]) == y) {
                add_one(bv_l, pos_to_add_l);
                pointers[x]++;
              } else {
                add_zero(bv_l, pos_to_add_l);
              }
            }
          }

#ifdef DEBUG
          cout << "Adding (..." << endl;
#endif // DEBUG
          add_one(bv_tree, pos_to_add);
#ifdef DEBUG
          cout << "Adding (..." << endl;
#endif // DEBUG
          add_one(bv_tree, pos_to_add);
#ifdef DEBUG
          cout << "Adding )..." << endl;
#endif // DEBUG
          add_zero(bv_tree, pos_to_add);
#ifdef DEBUG
          cout << "Adding )..." << endl;
#endif // DEBUG
          add_zero(bv_tree, pos_to_add);
        } else if(subm_size == k) {
#ifdef DEBUG
          cout << "Adding (..." << endl;
#endif // DEBUG
          add_one(bv_tree, pos_to_add);
#ifdef DEBUG
          cout << "Adding )..." << endl;
#endif // DEBUG
          add_zero(bv_tree, pos_to_add);
        } else {

          if(flag) {
#ifdef DEBUG
          cout << "Adding )..." << endl;
#endif // DEBUG
            add_zero(bv_tree, pos_to_add);
            continue;
          } else {
            recursion.push(make_tuple(subm_size, init_x, init_y, true));
          }

#ifdef DEBUG
          cout << "Adding (..." << endl;
#endif // DEBUG
          add_one(bv_tree, pos_to_add);
          
          if(has_at_least_one(adj_list, pointers, init_x, init_y, subm_size)) {
            recursion.push(make_tuple(subm_size / 2, init_x + subm_size / 2, init_y + subm_size / 2, false));
            recursion.push(make_tuple(subm_size / 2, init_x + subm_size / 2, init_y, false));
            recursion.push(make_tuple(subm_size / 2, init_x, init_y + subm_size / 2, false));
            recursion.push(make_tuple(subm_size / 2, init_x, init_y, false));
          } else {
#ifdef DEBUG
          cout << "Adding )..." << endl;
#endif // DEBUG
            add_zero(bv_tree, pos_to_add);
            recursion.pop();
          }
        }
      }
#ifdef DEBUG
      cout << "Adding ) of root ..." << endl;
#endif // DEBUG
      add_zero(bv_tree, pos_to_add);
      recursion.pop();

      last_bit_t = pos_to_add;
      last_bit_l = pos_to_add_l;

      tree = bit_vector(pos_to_add, 0);
      for(const auto& bit : bv_tree) tree[bit] = 1;

      l = bit_vector(pos_to_add_l, 0);
      for(const auto& bit : bv_l) l[bit] = 1;

      tree_support = bp_support_sada<>(&tree);
    }

    vector< pair< uint64_t, uint64_t > > get_pos_ones() {
      stack< tuple< uint8_t, uint64_t, uint64_t > > child_visit;

      uint64_t r, c;
      r = c = 0;
      child_visit.push({0, r, c});

      uint64_t to_read_l = 0;
      vector< pair< uint64_t, uint64_t > > ret;

      for(uint64_t i = 1; i < tree.size(); i++) {
#ifdef DEBUG
        cout << "Total bits " << tree.size() << endl;
        cout << "Reading " << i << " bit" << endl;
#endif // DEBUG
        if(tree[i]) {
#ifdef DEBUG
          cout << "Start of subtree" << endl;
#endif // DEBUG
          auto [vis, r_, c_] = child_visit.top();
          r = r_ + (vis / k) * (1 << (height_tree - child_visit.size()));
          c = c_ + (vis % k) * (1 << (height_tree - child_visit.size()));
          child_visit.push({0, r, c});
#ifdef DEBUG
          cout << "Level " << child_visit.size() << endl;
          cout << "Current row: " << r << " col: " << c << endl;
#endif // DEBUG
        } else {
#ifdef DEBUG
          cout << "End of subtree" << endl;
#endif // DEBUG
          if(child_visit.size() == height_tree + 1) {
#ifdef DEBUG
            cout << "Last level, read real values" << endl;
            cout << "L size " << l.size() << " ";
            cout << "Current bit " << to_read_l << endl;
#endif // DEBUG
            for(uint64_t j = 0; j < k * k; j++) {
              if(l[to_read_l]) {
                auto [vis, r_, c_] = child_visit.top();
                ret.push_back({r_ + j / k, c_ + j % k});
              }
              to_read_l++;
            }
#ifdef DEBUG
            cout << "Finishing reading real values" << endl;
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
      return ret;
    }

    void identical_trees() {
      csa_sada<> csa;
      string bp = "";
      for(uint64_t i = 0; i < tree.size(); i++) {
        bp += (tree[i] ? "(" : ")");
      }

      construct_im(csa, bp, 1);

      cout << csa << "\n";
      uint64_t amount_idem_subtree = 0;
      uint64_t amount_of_groups = 0;

      vector< uint64_t > pointer_or_not(tree.size(), -1);

      for(uint64_t pos_bp = 2; pos_bp < csa.size(); pos_bp++) {
        // only considering suffix starting with (
        if(tree[csa[pos_bp]]) {
          
          uint64_t curr_start_pos = csa[pos_bp];
          uint64_t curr_end_pos = tree_support.find_close(curr_start_pos);

          uint64_t prev_start_pos = csa[pos_bp - 1];
          uint64_t prev_end_pos = tree_support.find_close(prev_start_pos);

          // ignoring leaves
          if(curr_end_pos - curr_start_pos <= 3) continue;

          if(curr_end_pos - curr_start_pos == prev_end_pos - prev_start_pos) {
            if(curr_start_pos < prev_start_pos) {
              pointer_or_not[prev_start_pos] = (pointer_or_not[curr_start_pos] == -1 ? curr_start_pos : pointer_or_not[curr_start_pos]);
            } else {
              pointer_or_not[curr_start_pos] = (pointer_or_not[prev_start_pos] == -1 ? prev_start_pos : pointer_or_not[prev_start_pos]);
            }
            amount_idem_subtree++;
          } else {
#ifdef DEBUG
            cout << curr_start_pos << " " << curr_end_pos << "\n";

            for(uint64_t i = curr_start_pos; i < curr_end_pos + 1; i++) {
              cout << (tree[i] ? "(" : ")");
            }
            cout << endl;
#endif // DEBUG
            amount_of_groups++;
          }

        } else {
          // ( < ) = true
          break;
        }
      }
    }

    uint64_t size_in_bits() {
      return sizeof(uint64_t) * 5 +
             size_in_bytes(tree) * 8 +
             size_in_bytes(tree_support) * 8 +
             size_in_bytes(l) * 8;
    }

    friend ostream& operator<<(ostream& os, const k2tree_bp_sdsl<k> &k2tree) {
      cout << "Height Tree: " << k2tree.height_tree << endl;
      cout << "Tree: ";
      for(uint64_t i = 0; i < k2tree.tree.size(); i++) {
        cout << (k2tree.tree[i] ? "(" : ")");
      }
      cout << endl;
      cout << "L: ";
      for(uint64_t i = 0; i < k2tree.l.size(); i++) {
        if(i % 4 == 0 && !(i == 0)) cout << " ";
        cout << (k2tree.l[i] ? "1" : "0");
      }
      return os;
    }
};
#endif // !K2_TREE_BP_SDSL
