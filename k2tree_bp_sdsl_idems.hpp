#ifndef K2_TREE_BP_SDSL_IDEM
#define K2_TREE_BP_SDSL_IDEM
// std includes
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <utility>
#include <vector>

// local includes
#include "k2tree_bp_sdsl.hpp"
#include "libsais/include/libsais64.h"
#include "util.hpp"
//#include "Modificacion-S18/s18/head/s18_vector.hpp"

// sdsl includes
#include <sdsl/construct.hpp>
#include <sdsl/io.hpp>
#include <sdsl/rank_support_v5.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/sd_vector.hpp>
//#include <sdsl/rle_vector.hpp>
#include <sdsl/bp_support_sada.hpp>
#include <sdsl/util.hpp>
#include <sdsl/vectors.hpp>

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
  void clear() { e.clear(); }
};

// k2-tree
// parameters:
//   * k * k: amount of children per node
template< uint64_t k = 2, 
          class bit_vector_1 = bit_vector, class rank1_1 = rank_support_v5<>,
          class bit_vector_2 = bit_vector, class rank1_2 = rank_support_v5<>, class rank0_2 = rank_support_v5<0>,
                                           class select1_2 = select_support_mcl<>, class select0_2 = select_support_mcl<0> >
class k2tree_bp_sdsl_idems {
  private:
    //int_vector<> P;
    //vlc_vector<coder::fibonacci> P; // indo 5.2
    //dac_vector<> P; // indo 5.1
    dac_vector_dp<rrr_vector<127>> P;

    bit_vector_1 occ_PoL;
    rank1_1 rank1_occ_PoL;

    bit_vector_2 PoL;
    rank1_2 rank1_PoL;
    rank0_2 rank0_PoL;
    select1_2 select1_PoL;
    select0_2 select0_PoL;

    sd_vector<> real_tree;
    select_support_sd<> select_real_tree;

    uint64_t height_tree;

    bp_support_sada<> tree_support;
    bit_vector tree; // k2tree
    uint64_t last_bit_t; // universe

    rrr_vector<127> l; // real values
    uint64_t last_bit_l; // universe

    uint64_t msize;
    uint64_t rmsize;
    uint64_t m;

  public:
    uint64_t nodes() { return (tree_support.find_close(0) + 1) / 2; }

    k2tree_bp_sdsl_idems(k2tree_bp_sdsl<k> &k2tree) { 
      msize = k2tree.msize;
      rmsize = k2tree.rmsize;
      m = k2tree.m;

      height_tree = k2tree.height_tree;

      l = rrr_vector<127>(k2tree.l);

      //cout << "Building suffix array..." << endl;
      string bp = "";
      for(uint64_t i = 0; i < k2tree.tree.size(); i++) {
        bp += (k2tree.tree[i] ? "(" : ")");
      }

      uint8_t *text = reinterpret_cast<uint8_t*>(bp.data());
      int64_t *csa = new int64_t[bp.length()];
      int64_t *plcp = new int64_t[bp.length()];
      int64_t *lcp = new int64_t[bp.length()];

      if(libsais64(text, csa, bp.length(), 0, NULL) != 0) throw std::runtime_error("SA construction failed");

      if(libsais64_plcp(text, csa, plcp, bp.length()) != 0) throw std::runtime_error("PLCP array construction failed");

      if(libsais64_lcp(plcp, csa, lcp, bp.length()) != 0) throw std::runtime_error("LCP array construction failed");

      uint64_t amount_idem_subtree = 0;
      uint64_t amount_of_groups = 0;

      //cout << "Searching identical subtrees..." << endl;
      // remove later
      union_find idems_tree(k2tree.tree.size());

      for(uint64_t pos_bp = 2; pos_bp < bp.size(); pos_bp++) {
        // only considering suffix starting with (
        if(k2tree.tree[csa[pos_bp]]) {
          
          uint64_t curr_start_pos = csa[pos_bp];
          uint64_t curr_end_pos = k2tree.tree_support.find_close(curr_start_pos);

          uint64_t prev_start_pos = csa[pos_bp - 1];
          uint64_t prev_end_pos = k2tree.tree_support.find_close(prev_start_pos);

#ifdef DEBUG
          cout << "VS Subtree pos: " << prev_start_pos << " " << prev_end_pos << endl;
          cout << "   Subtree pos: " << curr_start_pos << " " << curr_end_pos << endl;
#endif
          // ignoring leaves
          if(curr_end_pos - curr_start_pos + 1 <= 4) continue;
          // ignoring small subtrees

          if(curr_end_pos - curr_start_pos == prev_end_pos - prev_start_pos) {
            bool flag = true;
            if(lcp[pos_bp] < curr_end_pos - curr_start_pos + 1) flag = false;

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
            cout << endl;
#endif // DEBUG
            amount_of_groups++;
          }

        } else {
          // ( < ) = true
          break;
        }
      }

      // clean, is useless
      free(csa);
      free(plcp);
      free(lcp);

      //cout << "Finding maximum head..." << endl;
      uint64_t maxi_repre = 0;
//      for(uint64_t nodes = 0; nodes < k2tree.tree.size(); nodes++) {
//        uint64_t repre = idems_tree.find_set(nodes);
//        if(k2tree.tree[nodes] == 1 && repre != nodes) {
//          maxi_repre = max(repre, maxi_repre);
//          nodes = k2tree.tree_support.find_close(nodes);
//        }
//      }
      vector< uint64_t > prefix_help(k2tree.tree.size(), 0);

      for(uint64_t bit = 0; bit < k2tree.tree.size(); bit++) {
        if(k2tree.tree[bit]) {

          uint64_t repre = idems_tree.find_set(bit);

          // if has an identical tree and is big enough
          if(repre != bit) {
            maxi_repre = (repre - prefix_help[repre - 1]);

            uint64_t next_bit = k2tree.tree_support.find_close(bit);
            for(uint64_t pfh = bit; pfh <= next_bit; pfh++) {
              prefix_help[pfh] = prefix_help[pfh - 1] + 1;
            }
            prefix_help[next_bit] -= 4;
            bit = next_bit;
          } else {
            if(bit > 0) prefix_help[bit] = prefix_help[bit - 1];
          }
        } else {
          prefix_help[bit] = prefix_help[bit - 1];
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
      uint64_t log2_w = ceil_log2(maxi_repre);

      //cout << "Replacing identical subtrees..." << endl;
      prefix_help.resize(k2tree.tree.size(), 0);

      for(uint64_t bit = 0; bit < k2tree.tree.size(); bit++) {
        if(k2tree.tree[bit]) {
          new_tree_bv.push_back(ref_bit++);

          uint64_t repre = idems_tree.find_set(bit);

          // if has an identical tree and is big enough
          if(repre != bit && k2tree.tree_support.find_close(repre) - repre + 1 > log2_w) {
            new_tree_bv.push_back(ref_bit++);
            pointer.push_back(repre - prefix_help[repre - 1]);
            ref_bit += 2;

            uint64_t next_bit = k2tree.tree_support.find_close(bit);
            for(uint64_t pfh = bit; pfh <= next_bit; pfh++) {
              prefix_help[pfh] = prefix_help[pfh - 1] + 1;
            }
            prefix_help[next_bit] -= 4;
            bit = next_bit;
          } else {
            if(bit > 0) prefix_help[bit] = prefix_help[bit - 1];
          }
        } else {
          ref_bit++;
          prefix_help[bit] = prefix_help[bit - 1];
        }
      }

      // clean, is useless
      prefix_help.clear();
      idems_tree.clear();


      tree = bit_vector(ref_bit, 0);
      for(const auto& bit : new_tree_bv) tree[bit] = 1;

      //cout << "Creating auxiliary bit vectors..." << endl;

      map< uint64_t, uint64_t > unique_pointer;
      uint64_t code = 0;

      for(uint64_t i = 0; i < pointer.size(); i++) {
        unique_pointer[pointer[i]] = 0;
      }

      for(auto itr = unique_pointer.begin(); itr != unique_pointer.end(); itr++) {
        itr->second = code++;
      }

      //cout << ceil_log2(code - 1) << " " << code << "\n";
      //cout << "idem subtrees: " << pointer.size() << "\n";

      int_vector<> aux(pointer.size());
      bit_vector bv_real_tree(tree.size(), 0);
      for(uint64_t i = 0; i < pointer.size(); i++) {
        //P[i] = unique_pointer[pointer[i]];
        aux[i] = unique_pointer[pointer[i]];
        bv_real_tree[pointer[i]] = 1;
      }

      real_tree = sd_vector(bv_real_tree);
      util::init_support(select_real_tree, &real_tree);

      // clean, is useless
      pointer.clear();
      // clean, is useless
      new_tree_bv.clear();

      vector< uint64_t > PoL_bv;
      vector< uint64_t > count_PoL;
      ref_bit = 0;

      //cout << "Creating PoL" << endl;

      for(uint64_t bit = 0, level = 0; bit < tree.size() - 3; bit++) {
        if(tree[bit]) level++;
        else level--;
        if(tree[bit] && tree[bit + 1] && !tree[bit + 2] && !tree[bit + 3]) {
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


      //util::bit_compress(aux);
      //P = vlc_vector<coder::fibonacci>(aux);
      //P = dac_vector<>(aux);
      P = dac_vector_dp<rrr_vector<127>>(aux);

      bit_vector bv_occ_PoL(count_PoL.back() + 1, 0);
      for(const auto& bit : count_PoL) bv_occ_PoL[bit] = 1;
      occ_PoL = bit_vector_1(bv_occ_PoL);
      util::init_support(rank1_occ_PoL, &occ_PoL);

      // clean, is useless
      count_PoL.clear();

      bit_vector bv_PoL(PoL_bv.back() + 1, 0);
      for(const auto& bit : PoL_bv) bv_PoL[bit] = 1;
      PoL = bit_vector_2(bv_PoL);

      util::init_support(rank1_PoL, &PoL);
      util::init_support(rank0_PoL, &PoL);
      util::init_support(select1_PoL, &PoL);
      util::init_support(select0_PoL, &PoL);

      // clean, is useless
      PoL_bv.clear();
      
      tree_support = bp_support_sada<>(&tree);

      //cout << "Finish building tree" << endl;
    }

    uint64_t access_PoL(uint64_t i) {
      uint64_t amount_ones = rank1_PoL(i + 1);
      if(amount_ones == 0) return 0;

      uint64_t pos_last_one = select1_PoL(amount_ones);

      if(pos_last_one < i) return 0;
      else return 1; 
    }

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

          if(P.size() > 0 && i < tree.size() - 4 && tree[i + 1] && !tree[i + 2] && !tree[i + 3]) {
            uint64_t read_PoL = (i >= occ_PoL.size() ? PoL.size() : rank1_occ_PoL(i));
#ifdef DEBUG
            cout << "Is a pointer or a leaf?" << endl;
            cout << "Current pos in PoL: " << read_PoL << endl;
            cout << "Size PoL: " << PoL.size() << endl;
#endif // DEBUG 
            if(read_PoL >= PoL.size() || access_PoL(read_PoL) == 0) continue;
#ifdef DEBUG
            cout << "Is a pointer" << endl;
#endif // DEBUG 

            uint64_t read_P = rank1_PoL(read_PoL);

#ifdef DEBUG
            cout << "To read from P: " << read_P << endl;
#endif // DEBUG 
            
            uint64_t where_to_move = select_real_tree(P[read_P] + 1);
            recover_pos.push({i, tree_support.find_close(where_to_move)});
            i = where_to_move;
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
            uint64_t value = l[to_read_l];
            for(uint64_t j = 0; j < k * k; j++) {
              if(value & (1 << j)) {
                auto [vis, r_, c_] = child_visit.top();
                ret.push_back({r_ + j / k, c_ + j % k});
              }
            }
            to_read_l++;
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

    uint64_t size_in_bits() {
      cout << "BITS" << endl;
      cout << "  Tree        : " << (size_in_bytes(tree)) * 8 << endl;
      cout << "  Tree Support: " << (size_in_bytes(tree_support)) * 8 << endl;
      cout << "  L           : " << (size_in_bytes(l)) * 8 << endl;
      cout << "  P           : " << (size_in_bytes(P)) * 8 << endl;
      cout << "  Real P      : " << (size_in_bytes(real_tree) + size_in_bytes(select_real_tree)) * 8 << endl;
      cout << "  occ_PoL     : " << (size_in_bytes(occ_PoL) + size_in_bytes(rank1_occ_PoL)) * 8 << endl;
      cout << "  PoL         : " << (size_in_bytes(PoL) + size_in_bytes(rank1_PoL) + size_in_bytes(rank0_PoL) + size_in_bytes(select1_PoL) + size_in_bytes(select0_PoL)) * 8 << endl;
      return sizeof(uint64_t) * 4 +
             size_in_bytes(tree) * 8 +
             size_in_bytes(tree_support) * 8 +
             size_in_bytes(l) * 8 + 
             size_in_bytes(P) * 8 +
             (size_in_bytes(real_tree) + size_in_bytes(select_real_tree)) * 8 +
             size_in_bytes(occ_PoL) * 8 + size_in_bytes(rank1_occ_PoL) * 8 +
             size_in_bytes(PoL) * 8 + size_in_bytes(rank1_PoL) * 8 + size_in_bytes(rank0_PoL) * 8 +
             size_in_bytes(select1_PoL) * 8 + size_in_bytes(select0_PoL) * 8;
    } 

    friend ostream& operator<<(ostream& os, const k2tree_bp_sdsl_idems< k, bit_vector_1, rank1_1, bit_vector_2, rank1_2, rank0_2, select1_2, select0_2 > &k2tree) {
      cout << "Height Tree: " << k2tree.height_tree << endl;
      cout << "Tree:      ";
      for(uint64_t i = 0; i < k2tree.tree.size(); i++) {
        cout << (k2tree.tree[i] ? "(" : ")");
      }
      cout << endl;
      cout << "L:         ";
      for(uint64_t i = 0; i < k2tree.l.size(); i++) {
        uint64_t value = k2tree.l[i];
        for(uint64_t j = 0; j < k * k; j++) {
          cout << (value & (1 << j) ? 1 : 0); 
        }
        cout << " ";
      }
      cout << endl;
      cout << "P:         ";
      for(uint64_t i = 0; i < k2tree.P.size(); i++) {
        cout << k2tree.P[i] << " ";
      }
      cout << endl;
      cout << "Real Tree: ";
      for(uint64_t i = 0; i < k2tree.real_tree.size(); i++) {
        cout << k2tree.real_tree[i];
      }
      cout << endl;
      cout << "PoL:       ";
      for(uint64_t i = 0; i < k2tree.PoL.size(); i++) {
        cout << k2tree.PoL[i];
      }
      cout << endl;
      cout << "OCC PoL:   ";
      for(uint64_t i = 0; i < k2tree.occ_PoL.size(); i++) {
        cout << k2tree.occ_PoL[i];
      }
      return os;
    }
};
#endif //K2_TREE_BP_SDSL_IDEM
