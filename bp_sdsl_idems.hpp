#ifndef TREE_BP_SDSL_IDEM
#define TREE_BP_SDSL_IDEM
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

// local includes
#include "libsais/include/libsais64.h"
#include "util.hpp"

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
  void clear() { e.clear(); }
};

// k2-tree
// parameters:
//   * k * k: amount of children per node
template< uint64_t k = 2 >
class bp_sdsl_idems {
  private:
    dac_vector_dp<rrr_vector<127>> P;

    sd_vector<> inside_idem;
    rank_support_sd<> rank1_inside_idem;

    // to replace inside_idem
    rrr_vector<127> one_or_two;
    sd_vector<> runs_0;

    sd_vector<> occ_PoL;
    rank_support_sd<> rank1_occ_PoL;

    sd_vector<> PoL;
    rank_support_sd<> rank1_PoL;
    rank_support_sd<0> rank0_PoL;
    select_support_sd<> select1_PoL;
    select_support_sd<0> select0_PoL;

    uint64_t height_tree;

    sd_vector<> real_tree;
    select_support_sd<> select_real_tree;

    bp_support_sada<> tree_support;
    bit_vector tree;
    uint64_t last_bit_t; // universe

  public:
    uint64_t nodes() { return (tree_support.find_close(0) + 1) / 2; }

    bp_sdsl_idems(bit_vector &tree) { 
      
      bp_support_sada<> tree_support(&tree);
      string bp = "";
      for(uint64_t i = 0; i < tree.size(); i++) {
        bp += (tree[i] ? "(" : ")");
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

      cout << "Searching identical subtrees..." << endl;
      // remove later
      union_find idems_tree(tree.size());

      for(uint64_t pos_bp = 2; pos_bp < bp.size(); pos_bp++) {
        // only considering suffix starting with (
        if(tree[csa[pos_bp]]) {
          
          uint64_t curr_start_pos = csa[pos_bp];
          uint64_t curr_end_pos = tree_support.find_close(curr_start_pos);

          uint64_t prev_start_pos = csa[pos_bp - 1];
          uint64_t prev_end_pos = tree_support.find_close(prev_start_pos);

#ifdef DEBUG
          cout << "VS Subtree pos: " << prev_start_pos << " " << prev_end_pos << endl;
          cout << "   Subtree pos: " << curr_start_pos << " " << curr_end_pos << endl;
#endif
          // ignoring leaves
          if(curr_end_pos - curr_start_pos + 1 <= 2) continue;
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

      cout << "Finding maximum head..." << endl;
      set< uint64_t > repres;
      vector< uint64_t > prefix_help(tree.size(), 0);

      for(uint64_t bit = 0; bit < tree.size(); bit++) {
        if(tree[bit]) {

          uint64_t repre = idems_tree.find_set(bit);

          // if has an identical tree and is big enough
          if(repre != bit) {
            repres.insert((repre - prefix_help[repre - 1]));

            uint64_t next_bit = tree_support.find_close(bit);
            for(uint64_t pfh = bit; pfh <= next_bit; pfh++) {
              prefix_help[pfh] = prefix_help[pfh - 1] + 1;
            }
            prefix_help[next_bit] -= 2;
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
      vector< uint64_t > inside_idem_bv;
      uint64_t ref_bit = 0;

      uint64_t amount_of_bits_removed = 0;
      uint64_t log2_w = ceil_log2(repres.size() - 1);
      cout << "Min pos " << *repres.begin() << endl;

      cout << "Replacing identical subtrees..." << endl;
      prefix_help.resize(tree.size(), 0);

      vector< uint64_t > pointers_pos;
      for(uint64_t bit = 0; bit < tree.size(); bit++) {
        if(tree[bit]) {
          new_tree_bv.push_back(ref_bit++);

          uint64_t repre = idems_tree.find_set(bit);

          // if has an identical tree and is big enough
          if(repre != bit && tree_support.find_close(repre) - repre + 1 > log2_w + 2) {
            inside_idem_bv.push_back(bit);
            pointers_pos.push_back(ref_bit - 1);
            //new_tree_bv.push_back(ref_bit++);
            pointer.push_back(repre - prefix_help[repre - 1]);
            ref_bit += 1;

            uint64_t next_bit = tree_support.find_close(bit);
            for(uint64_t pfh = bit; pfh <= next_bit; pfh++) {
              prefix_help[pfh] = prefix_help[pfh - 1] + 1;
            }
            prefix_help[next_bit] -= 2;
            bit = next_bit;
            inside_idem_bv.push_back(bit);
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

      vector< uint64_t > runs_0_pos;
      vector< uint64_t > one_or_two_pos;
      if(inside_idem_bv.size() > 0) {
        runs_0_pos.push_back(inside_idem_bv[0]);
      }


      uint64_t bit_12 = 0;
      for(uint64_t i = 1; i < inside_idem_bv.size(); i++) {
        if(inside_idem_bv[i] == inside_idem_bv[i - 1] + 1) {
          one_or_two_pos.push_back(bit_12);
        } else {
          runs_0_pos.push_back(runs_0_pos.back() + inside_idem_bv[i] - inside_idem_bv[i - 1]);
          bit_12++;
        }
      }

      auto aux_one_or_two = bit_vector(one_or_two_pos.back() + 1, 0);
      for(auto p : one_or_two_pos) aux_one_or_two[p] = 1;

      one_or_two = rrr_vector<127>(aux_one_or_two);

      runs_0 = sd_vector<>(runs_0_pos.begin(), runs_0_pos.end());

      bit_vector bv_aux(inside_idem_bv.back() + 1, 0);
      for(auto p : inside_idem_bv) bv_aux[p] = 1;

      inside_idem = sd_vector<>(bv_aux);
      sdsl::util::init_support(rank1_inside_idem, &inside_idem);

      this->tree = bit_vector(ref_bit, 0);
      for(const auto& bit : new_tree_bv) this->tree[bit] = 1;

      cout << "Creating auxiliary bit vectors..." << endl;

      map< uint64_t, uint64_t > unique_pointer;
      uint64_t code = 0;

      for(uint64_t i = 0; i < pointer.size(); i++) {
        unique_pointer[pointer[i]] = 0;
      }

      for(auto itr = unique_pointer.begin(); itr != unique_pointer.end(); itr++) {
        itr->second = code++;
      }

      cout << ceil_log2(code - 1) << " " << code << "\n";
      cout << "idem subtrees: " << pointer.size() << "\n";

      int_vector<> aux(pointer.size());
      bit_vector bv_real_tree(this->tree.size(), 0);
      for(uint64_t i = 0; i < pointer.size(); i++) {
        //P[i] = unique_pointer[pointer[i]];
        aux[i] = unique_pointer[pointer[i]];
        bv_real_tree[pointer[i]] = 1;
      }

      real_tree = sd_vector(bv_real_tree);
      sdsl::util::init_support(select_real_tree, &real_tree);

      // clean, is useless
      pointer.clear();
      // clean, is useless
      new_tree_bv.clear();

      vector< uint64_t > count_PoL;
      ref_bit = 0;

      cout << "Creating PoL" << endl;

      uint64_t order = 0;
      vector< uint64_t > PoL_bv;

      for(uint64_t bit = 0, level = 0; bit < this->tree.size() - 2; bit++) {
        if(this->tree[bit]) level++;
        else level--;
        if(this->tree[bit] && !this->tree[bit + 1]) {
#ifdef DEBUG
          cout << "Level: " << level << endl;
          cout << "Current bit: "  << bit << endl;
#endif
          count_PoL.push_back(bit);
          if(bit == pointers_pos[order]) {
            PoL_bv.push_back(ref_bit);
            order++;
          }
          ref_bit++;
        }
      }


      P = dac_vector_dp<rrr_vector<127>>(aux);

      bit_vector bv_occ_PoL((this->tree).size(), 0);
      for(const auto& bit : count_PoL) bv_occ_PoL[bit] = 1;
      occ_PoL = sd_vector<>(bv_occ_PoL);
      sdsl::util::init_support(rank1_occ_PoL, &occ_PoL);

      // clean, is useless
      count_PoL.clear();

      bit_vector bv_PoL(ref_bit, 0);
      for(const auto& bit : PoL_bv) bv_PoL[bit] = 1;
      PoL = sd_vector<>(bv_PoL);

      sdsl::util::init_support(rank1_PoL, &PoL);
      sdsl::util::init_support(rank0_PoL, &PoL);
      sdsl::util::init_support(select1_PoL, &PoL);
      sdsl::util::init_support(select0_PoL, &PoL);

      // clean, is useless
      PoL_bv.clear();
      
      this->tree_support = bp_support_sada<>(&(this->tree));

      cout << "Finish building tree" << endl;

    }


    uint64_t size_in_bits() {
      cout << "Tree           :" << size_in_bytes(tree) * 8 << endl;
      cout << "Tree support   :" <<  size_in_bytes(tree_support) * 8 << endl;
      cout << "P              :" << size_in_bytes(P) * 8 << endl;
      cout << "Real Tree      :" << size_in_bytes(real_tree) * 8 << endl;
      cout << "Inside idem    :" << size_in_bytes(inside_idem) * 8 + size_in_bytes(rank1_inside_idem) * 8 << endl;
      cout << "Inside runs_0  :" << size_in_bytes(runs_0) * 8 << endl; 
      cout << "Inside 1/2     :" << size_in_bytes(one_or_two) * 8 << endl;
      return sizeof(uint64_t) * 3 +
             size_in_bytes(tree) * 8 +
             size_in_bytes(tree_support) * 8 +
             size_in_bytes(P) * 8 +
             size_in_bytes(real_tree) * 8 +
             size_in_bytes(inside_idem) * 8 + size_in_bytes(rank1_inside_idem) * 8;
             //size_in_bytes(occ_PoL) * 8 + size_in_bytes(rank1_occ_PoL) * 8;
             //size_in_bytes(PoL) * 8 + size_in_bytes(rank1_PoL) * 8 + size_in_bytes(rank0_PoL) * 8 + size_in_bytes(select1_PoL) * 8 + size_in_bytes(select0_PoL) * 8;
    } 

    friend ostream& operator<<(ostream& os, const bp_sdsl_idems<k> &tree) {
      cout << "Tree:      ";
      for(uint64_t i = 0; i < tree.tree.size(); i++) {
        cout << (tree.tree[i] ? "(" : ")");
      }
      cout << endl;
      cout << "P:         ";
      for(uint64_t i = 0; i < tree.P.size(); i++) {
        cout << tree.P[i] << " ";
      }
      cout << endl;
      cout << "Real Tree: ";
      for(uint64_t i = 0; i < tree.real_tree.size(); i++) {
        cout << tree.real_tree[i];
      }
      cout << endl;
      cout << "PoL:       ";
      for(uint64_t i = 0; i < tree.PoL.size(); i++) {
        cout << tree.PoL[i];
      }
      cout << endl;
      cout << "OCC PoL:   ";
      for(uint64_t i = 0; i < tree.occ_PoL.size(); i++) {
        cout << tree.occ_PoL[i];
      }
      return os;
    }
};
#endif //K2_TREE_BP_SDSL_IDEM
