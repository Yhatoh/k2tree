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
#include "k2tree_bp_sdsl.hpp"

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
class k2tree_bp_sdsl_idems {
  private:
    int_vector<64> P;

    bit_vector PoL;
    // add rank/select

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
    k2tree_bp_sdsl_idems(k2tree_bp_sdsl<k> &k2tree) { 
      msize = k2tree.msize;
      rmsize = k2tree.rmsize;
      m = k2tree.m;

      height_tree = k2tree.height_tree;

      l = k2tree.l;

      csa_sada<> csa;
      string bp = "";
      for(uint64_t i = 0; i < k2tree.tree.size(); i++) {
        bp += (k2tree.tree[i] ? "(" : ")");
      }

      construct_im(csa, bp, 1);

      uint64_t amount_idem_subtree = 0;
      uint64_t amount_of_groups = 0;

      vector< uint64_t > pointer_or_not(k2tree.tree.size(), -1);

      for(uint64_t pos_bp = 2; pos_bp < csa.size(); pos_bp++) {
        // only considering suffix starting with (
        if(k2tree.tree[csa[pos_bp]]) {
          
          uint64_t curr_start_pos = csa[pos_bp];
          uint64_t curr_end_pos = k2tree.tree_support.find_close(curr_start_pos);

          uint64_t prev_start_pos = csa[pos_bp - 1];
          uint64_t prev_end_pos = k2tree.tree_support.find_close(prev_start_pos);

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

      vector< uint64_t > new_tree_bv;
      vector< uint64_t > pointer;
      uint64_t ref_bit = 0;

      for(uint64_t bit = 0; bit < k2tree.tree.size(); bit++) {
        if(k2tree.tree[bit]) {
          new_tree_bv.push_back(ref_bit++);

          if(pointer_or_not[bit] != -1) {
            new_tree_bv.push_back(ref_bit++);
            if(pointer_or_not[pointer_or_not[bit]] == -1) {
              pointer.push_back(pointer_or_not[bit]);
            } else {
              pointer.push_back(pointer_or_not[pointer_or_not[bit]]);
              pointer_or_not[bit] = pointer_or_not[pointer_or_not[bit]];
            }
            ref_bit += 2;

            bit = k2tree.tree_support.find_close(bit);
          }
        } else {
          ref_bit++;
        }
      }
      
      P = int_vector<64>(pointer.size());
      for(uint64_t i = 0; i < pointer.size(); i++) P[i] = pointer[i];

      tree = bit_vector(ref_bit, 0);
      for(const auto& bit : new_tree_bv) tree[bit] = 1;

      vector< uint64_t > PoL_bv;
      ref_bit = 0;
      for(uint64_t bit = 0, level = 0; bit < tree.size() - 3; bit++) {
        if(tree[bit]) level++;
        else level--;
        if(tree[bit] && tree[bit + 1] && !tree[bit + 2] && !tree[bit + 3]) {
          if(level != height_tree) {
            PoL_bv.push_back(ref_bit);
          }
          ref_bit++;
        }
      }

      PoL = bit_vector(ref_bit, 0);
      for(const auto& bit : PoL_bv) PoL[bit] = 1;

      tree_support = bp_support_sada<>(&tree);

    }

    vector< pair< uint64_t, uint64_t > > get_pos_ones() {
      return {};
    }

    void identical_trees() {
    }

    uint64_t size_in_bits() {
      return sizeof(uint64_t) * 5 +
             size_in_bytes(tree) * 8 +
             size_in_bytes(tree_support) * 8 +
             size_in_bytes(l) * 8;
    } 

    friend ostream& operator<<(ostream& os, const k2tree_bp_sdsl_idems<k> &k2tree) {
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
      cout << endl;
      cout << "P: ";
      for(uint64_t i = 0; i < k2tree.P.size(); i++) {
        cout << k2tree.P[i] << " ";
      }
      cout << endl;
      cout << "PoL: ";
      for(uint64_t i = 0; i < k2tree.PoL.size(); i++) {
        cout << k2tree.PoL[i];
      }
      cout << endl;
      return os;
    }
};
