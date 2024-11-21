#ifndef K2_TREE_BP_SDSL_IDEM
#define K2_TREE_BP_SDSL_IDEM
// std includes
#include <cassert>
#include <cstdint>
#include <iostream>
#include <sdsl/construct.hpp>
#include <sdsl/io.hpp>
#include <sdsl/lcp_support_sada.hpp>
#include <sdsl/rank_support_v5.hpp>
#include <utility>
#include <vector>

// sdsl includes
#include <sdsl/int_vector.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/bp_support_sada.hpp>
#include <sdsl/util.hpp>
#include <sdsl/lcp.hpp>

using namespace std;
using namespace sdsl;

struct union_find {
  vector< int64_t > e;
  union_find(int64_t n) { e.assign(n, -1); }
  int64_t find_set(int64_t x) { 
    return (e[x] < 0 ? x : e[x] = find_set(e[x]));
  }
  bool same_set (int64_t x, int64_t y) { return find_set(x) == find_set(y); }
  int64_t size (int64_t x) { return -e[find_set(x)]; }
  bool union_set (int64_t x, int64_t y) {
    x = find_set(x), y = find_set(y);
    if (x == y) return 0;
    if (x > y) swap(x, y);
    e[x] += e[y], e[y] = x;
    return 1;
  }
};

// k2-tree
// parameters:
//   * k * k: amount of children per node
template< uint64_t k = 2 >
class bp_sdsl_idems {
  private:
    int_vector<64> P;

    sd_vector<> occ_PoL;
    rank_support_sd<> rank1_occ_PoL;

    sd_vector<> PoL;
    rank_support_sd<> rank1_PoL;
    rank_support_sd<0> rank0_PoL;
    select_support_sd<> select1_PoL;
    select_support_sd<0> select0_PoL;

    uint64_t height_tree;

    bp_support_sada<> tree_support;
    bit_vector tree; // k2tree
    uint64_t last_bit_t; // universe

  public:
    uint64_t nodes() { return (tree_support.find_close(0) + 1) / 2; }

    bp_sdsl_idems(bit_vector &bp) { 
      
      csa_sada<> csa;
      string string_bp = "";
      for(uint64_t i = 0; i < bp.size(); i++) {
        string_bp += (bp[i] ? "(" : ")");
      }
      

#ifdef DEBUG
      cout << "Building suffix array" << endl;
#endif

      construct_im(csa, string_bp, 1);

#ifdef DEBUG
      cout << "Building tree support" << endl;
#endif
      bp_support_sada<> tree_support_aux(&bp);

#ifdef DEBUG
      cout << "Starting searching identical subtrees" << endl;
#endif
      uint64_t amount_idem_subtree = 0;
      uint64_t amount_of_groups = 0;

      // remove later
      union_find idems_tree(bp.size());


      for(uint64_t pos_bp = 2; pos_bp < csa.size(); pos_bp++) {
        // only considering suffix starting with (
        if(bp[csa[pos_bp]]) {
          
          uint64_t curr_start_pos = csa[pos_bp];
          uint64_t curr_end_pos = tree_support_aux.find_close(curr_start_pos);

          uint64_t prev_start_pos = csa[pos_bp - 1];
          uint64_t prev_end_pos = tree_support_aux.find_close(prev_start_pos);

#ifdef DEBUG
          cout << "VS Subtree pos: " << prev_start_pos << " " << prev_end_pos << endl;
          cout << "   Subtree pos: " << curr_start_pos << " " << curr_end_pos << endl;
#endif
          // ignoring leaves
          if(curr_end_pos - curr_start_pos <= 2) continue;

          if(curr_end_pos - curr_start_pos == prev_end_pos - prev_start_pos) {
            bool flag = true;
            for(uint64_t i = 0; i < curr_end_pos - curr_start_pos + 1; i++) {
              if(bp[curr_start_pos + i] != bp[prev_start_pos + i]) {
                flag = false;
                break;
              }
            }

            if(flag) {
              idems_tree.union_set(curr_start_pos, prev_start_pos);
#ifdef DEBUG
              cout << "== Subtree pos: " << prev_start_pos << " " << prev_end_pos << endl;
              cout << "   Subtree pos: " << curr_start_pos << " " << curr_end_pos << endl;
              cout << "   Head: " << idems_tree.find_set(curr_start_pos) << endl;
#endif
              amount_idem_subtree++;
            } else {
#ifdef DEBUG
              cout << "!= Subtree pos: " << prev_start_pos << " " << prev_end_pos << endl;
              cout << "   Subtree pos: " << curr_start_pos << " " << curr_end_pos << endl;
#endif
              amount_of_groups++;
            }
          } else {
#ifdef DEBUG
            cout << "!= Subtree pos: " << prev_start_pos << " " << prev_end_pos << endl;
            cout << "   Subtree pos: " << curr_start_pos << " " << curr_end_pos << endl;
            cout << "Subtree Start: " << curr_start_pos << " End: " << curr_end_pos << endl;
            //for(uint64_t i = curr_start_pos; i < curr_end_pos + 1; i++) {
            //  cout << (k2tree.tree[i] ? "(" : ")");
            //}
            cout << endl;
#endif // DEBUG
            amount_of_groups++;
          }

        } else {
          // ( < ) = true
          break;
        }
      }

#ifdef DEBUG
      cout << amount_idem_subtree << endl;
      cout << amount_of_groups << endl;
#endif // DEBUG
      
      vector< uint64_t > new_tree_bv;
      vector< uint64_t > pointer;
      uint64_t ref_bit = 0;

      uint64_t amount_of_bits_removed = 0;

      vector< uint64_t > prefix_help(bp.size(), 0);

      for(uint64_t bit = 0; bit < bp.size(); bit++) {
        if(bp[bit]) {
          new_tree_bv.push_back(ref_bit++);

          uint64_t repre = idems_tree.find_set(bit);

          if(repre != bit) {
            pointer.push_back(repre - prefix_help[repre - 1]);
            ref_bit += 1;

            uint64_t next_bit = tree_support_aux.find_close(bit);
            for(uint64_t pfh = bit; pfh <= next_bit; pfh++) {
              prefix_help[pfh] = prefix_help[pfh - 1] + 1;
            }
            prefix_help[next_bit] -= 2;
            bit = next_bit;
          } else {
            if(bit > 0) prefix_help[bit] = prefix_help[bit - 1];
          }
        } else {
          ref_bit++;
          prefix_help[bit] = prefix_help[bit - 1];
        }

      }
      
      P = int_vector<64>(pointer.size());
      for(uint64_t i = 0; i < pointer.size(); i++) P[i] = pointer[i];

      tree = bit_vector(ref_bit, 0);
      for(const auto& bit : new_tree_bv) tree[bit] = 1;

      vector< uint64_t > PoL_bv;
      vector< uint64_t > count_PoL;
      ref_bit = 0;
#ifdef DEBUG
      cout << "Creating PoL" << endl;
#endif
      for(uint64_t bit = 0, level = 0; bit < tree.size() - 3; bit++) {
        if(tree[bit]) level++;
        else level--;
        if(tree[bit] && !tree[bit + 1]) {
#ifdef DEBUG
          cout << "Level: " << level << endl;
          cout << "Current bit: "  << bit << endl;
#endif
          count_PoL.push_back(bit);
          if(level != height_tree) {
            PoL_bv.push_back(ref_bit);
          }
          ref_bit++;
        }
      }

      cout << count_PoL.size() << "\n";

      occ_PoL = sd_vector<>(count_PoL.begin(), count_PoL.end());
      util::init_support(rank1_occ_PoL, &occ_PoL);

      PoL = sd_vector<>(PoL_bv.begin(), PoL_bv.end());
      util::init_support(rank1_PoL, &PoL);
      util::init_support(rank0_PoL, &PoL);
      util::init_support(select1_PoL, &PoL);
      util::init_support(select0_PoL, &PoL);
      
      tree_support = bp_support_sada<>(&tree);

    }

    /*
    vector< pair< uint64_t, uint64_t > > get_pos_ones() {
      stack< tuple< uint8_t, uint64_t, uint64_t > > child_visit;
      stack< pair< uint64_t, uint64_t > > recover_pos;

      uint64_t r, c;
      r = c = 0;
      child_visit.push({0, r, c});

      uint64_t to_read_l = 0;

      vector< pair< uint64_t, uint64_t > > ret;

      for(uint64_t i = 1; i < tree.size(); i++) {
#ifdef DEBUG
        cout << "-------------------" << endl;
        cout << "Total bits " << tree.size() << endl;
        cout << "Reading " << i << " bit" << endl;
#endif // DEBUG
        if(tree[i]) {
#ifdef DEBUG
          cout << "Start of subtree" << endl;
#endif // DEBUG
          
          auto [vis, r_, c_] = child_visit.top();
#ifdef DEBUG
          cout << "Previous tree: " << (uint64_t) vis << " " << r_ << " " << c_ << endl;
#endif // DEBUG
          r = r_ + (vis / k) * (1 << (height_tree - child_visit.size()));
          c = c_ + (vis % k) * (1 << (height_tree - child_visit.size()));
          child_visit.push({0, r, c});

          if(i < tree.size() - 4 && tree[i + 1] && !tree[i + 2] && !tree[i + 3]) {
            uint64_t read_PoL = (i >= occ_PoL.size() ? PoL.size() : rank1_occ_PoL(i));
#ifdef DEBUG
            cout << "Is a pointer or a leaf?" << endl;
            cout << "Current pos in PoL: " << read_PoL << endl;
            cout << "Size PoL: " << PoL.size() << endl;
#endif // DEBUG 
            if(read_PoL >= PoL.size() || PoL[read_PoL] == 0) continue;
#ifdef DEBUG
            cout << "Is a pointer" << endl;
#endif // DEBUG 

            uint64_t read_P = rank1_PoL(read_PoL);

#ifdef DEBUG
            cout << "To read from P: " << read_P << endl;
#endif // DEBUG 
            
            recover_pos.push({i, tree_support.find_close(P[read_P])});
            i = P[read_P];
#ifdef DEBUG
            cout << "Moving to pos: " << i << endl;
#endif // DEBUG 
          }
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
          if(!recover_pos.empty() && recover_pos.top().second <= i) {
#ifdef DEBUG
            cout << "Finish reading copy, returning to original tree at pos: " << recover_pos.top().first << endl;
#endif // DEBUG
            i = recover_pos.top().first + 3;
            recover_pos.pop();
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
*/

    uint64_t size_in_bits() {
      return sizeof(uint64_t) * 5 +
             size_in_bytes(tree) * 8 +
             size_in_bytes(tree_support) * 8 +
             size_in_bytes(P) * 8 +
             size_in_bytes(occ_PoL) * 8 + size_in_bytes(rank1_occ_PoL) * 8;
             //size_in_bytes(PoL) * 8 + size_in_bytes(rank1_PoL) * 8 + size_in_bytes(rank0_PoL) * 8 + size_in_bytes(select1_PoL) * 8 + size_in_bytes(select0_PoL) * 8;
    } 

    friend ostream& operator<<(ostream& os, const bp_sdsl_idems<k> &tree) {
      cout << "Height Tree: " << tree.height_tree << endl;
      cout << "Tree: ";
      for(uint64_t i = 0; i < tree.tree.size(); i++) {
        cout << (tree.tree[i] ? "(" : ")");
      }
      cout << endl;
      cout << "P: ";
      for(uint64_t i = 0; i < tree.P.size(); i++) {
        cout << tree.P[i] << " ";
      }
      cout << endl;
      cout << "PoL: ";
      for(uint64_t i = 0; i < tree.PoL.size(); i++) {
        cout << tree.PoL[i];
      }
      cout << endl;
      cout << "OCC PoL: ";
      for(uint64_t i = 0; i < tree.occ_PoL.size(); i++) {
        cout << tree.occ_PoL[i];
      }
      return os;
    }
};
#endif //K2_TREE_BP_SDSL_IDEM
