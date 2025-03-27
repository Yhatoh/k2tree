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
#include "plaintree.hpp"
#include "util.hpp"
//#include "Modificacion-S18/s18/head/s18_vector.hpp"

// sdsl includes
#include <sdsl/construct.hpp>
#include <sdsl/io.hpp>
#include <sdsl/rank_support_v5.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/bp_support_sada.hpp>
#include <sdsl/util.hpp>
#include <sdsl/vectors.hpp>

using namespace std;
using namespace sdsl;

#define dbg(var) cout << #var << " = " << var << endl;

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
          class bv_leaves = bit_vector,
          class bit_vector_1 = bit_vector, class rank1_1 = rank_support_v5<>,
          class bit_vector_2 = bit_vector, class rank1_2 = rank_support_v5<>, class rank0_2 = rank_support_v5<0>,
                                           class select1_2 = select_support_mcl<>, class select0_2 = select_support_mcl<0>, uint64_t b_size = 1280 >
class k2tree_bp_sdsl_idems {
  private:
    //int_vector<> P;
    //vlc_vector<coder::fibonacci> P; // indo 5.2
    //dac_vector<> P; // indo 5.1
    dac_vector_dp<rrr_vector<127>> P;

    //bit_vector_1 occ_PoL;
    //rank1_1 rank1_occ_PoL;

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

    bv_leaves l; // real values
    uint64_t last_bit_l; // universe

    uint64_t msize;
    uint64_t rmsize;
    uint64_t m;

    uint64_t maximal_subtrees; // just for information purposes

    int_vector<> blocks;

    void add_one(vector< uint64_t > &bv, uint64_t &pos_to_add) {
      bv.push_back(pos_to_add++);
    }

    void add_zero(vector< uint64_t > &bv, uint64_t &pos_to_add) {
      pos_to_add++;
    }

    inline uint64_t count(uint64_t num) {
      uint64_t x = num;
      uint64_t y = num >> 1;
      uint64_t hi = x & y;
      uint64_t lo = ~ (x | y);

      uint64_t bits = (hi & (lo >> 2)) & 2305843009213693951;
      return __builtin_popcountll(bits);
    }

    inline uint64_t rank_occ_PoL(uint64_t i) {
      uint64_t b = i / b_size;
      uint64_t ret = (b > 0 ? blocks[b - 1] : 0);

      uint64_t bit = b * b_size;
      for(; bit + 64 < i; bit += 64) {
        uint64_t extra = 0;
        uint64_t read = tree.get_int(bit, 64);
        if(bit + 64 < i) {
          uint64_t len = (6 > i - (bit + 61) ? i - (bit + 61) : 6);
          extra = count(tree.get_int(bit + 61, len) | (((uint64_t) -1) << len));
        }
        ret += count(read) + extra;
      }
      if(i > bit) {
        if(i - bit == 64) ret += count(tree.get_int(bit, 64));
        else ret += count(tree.get_int(bit, i - bit) | (((uint64_t) -1) << (i - bit)));
      }
      return ret;
    }

  public:
    uint64_t height() { return height_tree; } 
    uint64_t size() { return m; }
    uint64_t size_matrix() { return rmsize; }
    uint64_t size_comp_subtrees() { return P.size(); }
    uint64_t size_maximal_subtrees() { return maximal_subtrees; }
    uint64_t nodes() { return (tree_support.find_close(0) + 1) / 2; }

    k2tree_bp_sdsl_idems() {}

    k2tree_bp_sdsl_idems(k2tree_bp_sdsl<k, bv_leaves> &k2tree) {
      k2tree.tree_support = bp_support_sada<>(&k2tree.tree);
      msize = k2tree.msize;
      rmsize = k2tree.rmsize;
      m = k2tree.m;

      height_tree = k2tree.height_tree;

      l = k2tree.l;

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
      set< uint64_t > maxi_repre;
      vector< uint64_t > prefix_help(k2tree.tree.size(), 0);


      for(uint64_t bit = 0; bit < k2tree.tree.size(); bit++) {
        if(k2tree.tree[bit]) {

          uint64_t repre = idems_tree.find_set(bit);

          // if has an identical tree and is big enough
          if(repre != bit) {
            //maxi_repre = (repre - prefix_help[repre - 1]);
            maxi_repre.insert(repre - prefix_help[repre - 1]);

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
      uint64_t log2_w = ceil_log2(maxi_repre.size() - 1);

      //cout << "Replacing identical subtrees..." << endl;
      prefix_help.resize(k2tree.tree.size(), 0);

      for(uint64_t bit = 0; bit < k2tree.tree.size(); bit++) {
        if(k2tree.tree[bit]) {
          new_tree_bv.push_back(ref_bit++);

          uint64_t repre = idems_tree.find_set(bit);

          // if has an identical tree and is big enough
          if(repre != bit && k2tree.tree_support.find_close(repre) - repre + 1 > log2_w + 4) {
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

      maximal_subtrees = code;

      int_vector<> aux(pointer.size());
      bit_vector bv_real_tree(tree.size(), 0);
      for(uint64_t i = 0; i < pointer.size(); i++) {
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


      //P = vlc_vector<coder::fibonacci>(aux);
      //P = dac_vector<>(aux);
      P = dac_vector_dp<rrr_vector<127>>(aux);

      bit_vector bv_occ_PoL(tree.size(), 0);
      for(const auto& bit : count_PoL) bv_occ_PoL[bit] = 1;
      occ_PoL = bit_vector_1(bv_occ_PoL);
      util::init_support(rank1_occ_PoL, &occ_PoL);

      // clean, is useless
      count_PoL.clear();

      bit_vector bv_PoL(ref_bit, 0);
      for(const auto& bit : PoL_bv) bv_PoL[bit] = 1;
      PoL = bit_vector_2(bv_PoL);

      util::init_support(rank1_PoL, &PoL);
      util::init_support(rank0_PoL, &PoL);
      util::init_support(select1_PoL, &PoL);
      util::init_support(select0_PoL, &PoL);

      // clean, is useless
      PoL_bv.clear();

      tree_support = bp_support_sada<>(&tree);

       // count (()) prefix sum
      uint64_t sum = 0;
      vector< uint64_t > prefix_sum;
      for(uint64_t bit = 1; bit < tree.size() - 3; bit++) {
        if(bit % b_size == 0) {
          prefix_sum.push_back(sum);
          //sum = 0;
        }
        if(tree[bit] && tree[bit + 1] && !tree[bit + 2] && !tree[bit + 3]) {
          sum += 1;
        }
      }

      prefix_sum.push_back(sum);
      blocks = int_vector<>(prefix_sum.size(), -1);
      for(uint64_t block = 0; block < prefix_sum.size(); block++) {
        blocks[block] = prefix_sum[block];
      }
      util::bit_compress(blocks);
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
            //uint64_t read_PoL = (i >= occ_PoL.size() ? PoL.size() : rank1_occ_PoL(i));
            uint64_t read_PoL = rank_occ_PoL(i);
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

    /*
    void binsum(const k2tree_bp_sdsl_idems< k, bv_leaves, bit_vector_1, rank1_1, bit_vector_2, rank1_2, rank0_2, select1_2, select0_2 >& B, plain_tree &C) {
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
        C.l = bit_vector(l.size(), 0);
        for(uint64_t i = 0; i < l.size(); i++) C.l[i] = l[i];
        C.height_tree = height_tree;
        C.msize = msize;
        C.rmsize = rmsize;
        C.m = m;
        return;
      }

      if(tree.size() == 2) {
        C.tree = B.tree;
        C.l = bit_vector(B.l.size(), 0);
        for(uint64_t i = 0; i < B.l.size(); i++) C.l[i] = B.l[i];
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

      stack< pair< uint64_t, uint64_t > > recover_pos_A;
      stack< pair< uint64_t, uint64_t > > recover_pos_B;
#ifdef DEBUG
      cout << "Starting algorithm" << endl;
#endif

      int64_t curr_depth = 0;
      while(pa < tree.size() && pb < B.tree.size()) {
        // see if there is a pointer A
        if(P.size() > 0 && pa < tree.size() - 4 && tree[pa] && tree[pa + 1] && !tree[pa + 2] && !tree[pa + 3]) {
          uint64_t read_PoL = (pa >= occ_PoL.size() ? PoL.size() : rank1_occ_PoL(pa));
          if(!(read_PoL >= PoL.size() || access_PoL(read_PoL) == 0)) {
            uint64_t read_P = rank1_PoL(read_PoL);

            uint64_t where_to_move = select_real_tree(P[read_P] + 1);
            recover_pos_A.push({pa, tree_support.find_close(where_to_move)});
            pa = where_to_move;
          }
        }
        // see if there is a pointer B
        if(B.P.size() > 0 && pb < B.tree.size() - 4 && B.tree[pb] && B.tree[pb + 1] && !B.tree[pb + 2] && !B.tree[pb + 3]) {
            uint64_t read_PoL = (pb >= B.occ_PoL.size() ? B.PoL.size() : B.rank1_occ_PoL(pb));
            if(!(read_PoL >= B.PoL.size() || B.PoL[read_PoL] == 0)) {
              uint64_t read_P = B.rank1_PoL(read_PoL);

              uint64_t where_to_move = B.select_real_tree(B.P[read_P] + 1);
              recover_pos_B.push({pb, B.tree_support.find_close(where_to_move)});
              pb = where_to_move;
            }
        }

        // check if you need to comeback A
        if(!recover_pos_A.empty() && recover_pos_A.top().second <= pa) {
          pa = recover_pos_A.top().first + 3;
          recover_pos_A.pop();
        }

        // check if you need to comeback B
        if(!recover_pos_B.empty() && recover_pos_B.top().second <= pb) {
          pb = recover_pos_B.top().first + 3;
          recover_pos_B.pop();
        }

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
               (pLb < B.l.size() && B.l[pLb]))
              add_one(bits_L, curr_bit_L);
            else
             add_zero(bits_L, curr_bit_L);
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
            // see if there is a pointer A
            if(P.size() > 0 && pa < tree.size() - 4 && tree[pa] && tree[pa + 1] && !tree[pa + 2] && !tree[pa + 3]) {
              uint64_t read_PoL = (pa >= occ_PoL.size() ? PoL.size() : rank1_occ_PoL(pa));
              if(!(read_PoL >= PoL.size() || access_PoL(read_PoL) == 0)) {
                uint64_t read_P = rank1_PoL(read_PoL);

                uint64_t where_to_move = select_real_tree(P[read_P] + 1);
                recover_pos_A.push({pa, tree_support.find_close(where_to_move)});
                pa = where_to_move;
              }
            }

            // check if you need to comeback A
            if(!recover_pos_A.empty() && recover_pos_A.top().second <= pa) {
              pa = recover_pos_A.top().first + 3;
              recover_pos_A.pop();
            }

            if(tree[pa] && curr_depth < height_tree) {
              add_one(bits_tree, curr_bit_tree);
              counter++;
              curr_depth++;
            } else if(tree[pa]) {
              add_one(bits_tree, curr_bit_tree);
              for(uint64_t i = 0; i < 4; i++) {
                if(pLa < l.size() && l[pLa]) add_one(bits_L, curr_bit_L);
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
            // see if there is a pointer B
            if(B.P.size() > 0 && pb < B.tree.size() - 4 && B.tree[pb] && B.tree[pb + 1] && !B.tree[pb + 2] && !B.tree[pb + 3]) {
              uint64_t read_PoL = (pb >= B.occ_PoL.size() ? B.PoL.size() : B.rank1_occ_PoL(pb));
              if(!(read_PoL >= B.PoL.size() || B.PoL[read_PoL] == 0)) {
                uint64_t read_P = B.rank1_PoL(read_PoL);

                uint64_t where_to_move = B.select_real_tree(B.P[read_P] + 1);
                recover_pos_B.push({pb, B.tree_support.find_close(where_to_move)});
                pb = where_to_move;
              }
            }

            // check if you need to comeback B
            if(!recover_pos_B.empty() && recover_pos_B.top().second <= pb) {
              pb = recover_pos_B.top().first + 3;
              recover_pos_B.pop();
            }

            if(B.tree[pb] && curr_depth < B.height_tree) {
              add_one(bits_tree, curr_bit_tree);
              counter++;
              curr_depth++;
            } else if(B.tree[pb]) {
              add_one(bits_tree, curr_bit_tree);
              for(uint64_t i = 0; i < 4; i++) {
                if(pLb < B.l.size() && B.l[pLb++]) add_one(bits_L, curr_bit_L);
                else add_zero(bits_L, curr_bit_L);
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
*/

    void prefix_sum_skipped_values(vector< uint64_t > &pre_skips) {
      uint64_t leaves = 0;
      uint64_t id = 0;

      stack< uint8_t > child_visit;
      stack< tuple< uint64_t, uint64_t, uint64_t, uint64_t > > recover_pos;
      vector< uint8_t > flags(pre_skips.size(), 0);

      child_visit.push(0);

      for(uint64_t i = 1; i < tree.size(); i++) {
#ifdef DEBUG
        cout << "-------------------" << endl;
        cout << "Total bits " << tree.size() << endl;
        cout << "Current Level " << child_visit.size() << endl;
        cout << "Pointers? " << recover_pos.size() << endl;
        cout << "Reading " << i << " bit" << endl;
#endif // DEBUG
        if(tree[i]) {
#ifdef DEBUG
          cout << "Start of subtree" << endl;
#endif // DEBUG

          child_visit.push(0);

          if(P.size() > 0 && i < tree.size() - 4 && tree[i + 1] && !tree[i + 2] && !tree[i + 3]) {
            uint64_t read_PoL = rank_occ_PoL(i);
            //uint64_t read_PoL = (i >= occ_PoL.size() ? PoL.size() : rank1_occ_PoL(i));
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
            recover_pos.push({i, tree_support.find_close(where_to_move), read_P, 0});
            i = where_to_move;
#ifdef DEBUG
            cout << "Moving to pos: " << i << endl;
#endif // DEBUG
          }
#ifdef DEBUG
          cout << "Level " << child_visit.size() << endl;
#endif // DEBUG
        } else {
#ifdef DEBUG
          cout << "End of subtree" << endl;
#endif // DEBUG
          if(child_visit.size() == height_tree + 1) {
#ifdef DEBUG
            cout << "Last level sum" << endl;
#endif // DEBUG
            if(!recover_pos.empty()) {
              auto [index, end, _id, leaves] = recover_pos.top();
              leaves += 4;
              recover_pos.pop();
              recover_pos.push({index, end, _id, leaves});
            }
#ifdef DEBUG
            cout << "Finishing reading real values" << endl;
#endif // DEBUG
          }
          if(!recover_pos.empty() && get<1>(recover_pos.top()) <= i) {
            auto [index, end, _id, leaves] = recover_pos.top();
            if(flags[_id] == 0) {
              pre_skips[_id] = leaves;
              flags[_id] = 1;
            }
            i = index + 3;

            recover_pos.pop();
            if(!recover_pos.empty()) {
              auto [s1, s2, s3, s4] = recover_pos.top();
              recover_pos.pop();
              recover_pos.push({s1, s2, s3, leaves + s4});
            }

          }
          child_visit.pop();
          // means we finish to read the complete tree
          if(child_visit.size() != 0) {
            auto vis = child_visit.top();
            child_visit.pop();
            child_visit.push(vis + 1);
          }
        }
      }

      for(uint64_t i = 1; i < pre_skips.size(); i++) {
        pre_skips[i] = pre_skips[i] + pre_skips[i - 1];
      }
    }

    void mul(k2tree_bp_sdsl_idems< k, bv_leaves, bit_vector_1, rank1_1, bit_vector_2, rank1_2, rank0_2, select1_2, select0_2 >& B, plain_tree &C) {
      vector< uint64_t > pre_skips_A(rank1_PoL(PoL.size()), 0);
      prefix_sum_skipped_values(pre_skips_A);

      vector< uint64_t > pre_skips_B(B.rank1_PoL(B.PoL.size()), 0);
      B.prefix_sum_skipped_values(pre_skips_B);

      sdsl::int_vector<4> A_L_S(l.size() / 4, 0);
      sdsl::int_vector<4> B_L_S(B.l.size() / 4, 0);

      uint64_t A_tree, B_tree;
      A_tree = B_tree = 0;
      uint64_t A_L, B_L;
      A_L = B_L = 0;
      mul(A_tree, A_L, A_L_S, pre_skips_A, 0, 0, B,
          B_tree, B_L, B_L_S, pre_skips_B, 0, 0, C, height_tree);
    }

    void mul(uint64_t &A_tree,
             uint64_t &A_L, sdsl::int_vector<4> &A_L_S,
             vector< uint64_t > &pre_skips_A, uint64_t A_lvs_sk, bool A_flag,
             k2tree_bp_sdsl_idems< k, bv_leaves, bit_vector_1, rank1_1,
                                         bit_vector_2, rank1_2, rank0_2, select1_2, select0_2 >& B,
             uint64_t &B_tree,
             uint64_t &B_L, sdsl::int_vector<4> &B_L_S,
             vector< uint64_t > &pre_skips_B, uint64_t B_lvs_sk, bool B_flag,
             plain_tree &C,
             uint64_t curr_h) {
#ifdef DEBUG
      cout << "Current Height: " << curr_h << endl;
      cout << "  Start Matrix A: " << A_tree << endl;
      cout << "  Pos in L     A: " << A_L << endl;
      cout << "  Start Matrix B: " << B_tree << endl;
      cout << "  Pos in L     B: " << B_L << endl;
#endif
      // submatrix A or B full of 0's
      bool A_f0 = (!tree[A_tree + 1]);
      bool B_f0 = (!B.tree[B_tree + 1]);
      if(A_f0 && B_f0) {
#ifdef DEBUG
        cout << "Full of 0's" << endl;
#endif
        C.reserve(2, 0);
        C.tree.push_back(1);
        C.tree.push_back(0);
        C.height_tree = curr_h;
        C.m = m;
        C.msize = msize;
        C.rmsize = rmsize;
        A_tree++;
        B_tree++;
        return;
      } else if(A_f0) {
#ifdef DEBUG
        cout << "A Full of 0's" << endl;
#endif
        C.reserve(2, 0);
        C.tree.push_back(1);
        C.tree.push_back(0);
        C.height_tree = curr_h;
        C.m = m;
        C.msize = msize;
        C.rmsize = rmsize;
        A_tree++;
        if(B_flag) return;
        B_tree = B.tree_support.find_close(B_tree);
        //uint64_t rank_p_B = B.rank1_occ_PoL(B_tree);
        uint64_t rank_p_B = B.rank_occ_PoL(B_tree);
        uint64_t pointer = B.rank1_PoL(rank_p_B);
        B_L = B_lvs_sk + B.rank0_PoL(rank_p_B) * 4 + (pointer > 0 ? pre_skips_B[pointer - 1] : 0);
        return;
      } else if(B_f0) {
#ifdef DEBUG
        cout << "B Full of 0's" << endl;
#endif
        C.reserve(2, 0);
        C.tree.push_back(1);
        C.tree.push_back(0);
        C.height_tree = curr_h;
        C.m = m;
        C.msize = msize;
        C.rmsize = rmsize;
        B_tree++;
        if(A_flag) return;

        A_tree = tree_support.find_close(A_tree);
        //uint64_t rank_p_A = rank1_occ_PoL(A_tree);
        uint64_t rank_p_A = rank_occ_PoL(A_tree);
        uint64_t pointer = rank1_PoL(rank_p_A);
        A_L = A_lvs_sk + rank0_PoL(rank_p_A) * 4 + (pointer > 0 ? pre_skips_A[pointer - 1] : 0);
        return;
      }

      // base case, leaf
      if(curr_h == 1) {
#ifdef DEBUG
        cout << "Leaf!" << endl;
#endif

        uint8_t aux_l = minimat_mul((A_L_S[A_L >> 2] ? A_L_S[A_L >> 2] : A_L_S[A_L >> 2] = l.get_int(A_L, 4)),
                                    (B_L_S[B_L >> 2] ? B_L_S[B_L >> 2] : B_L_S[B_L >> 2] = B.l.get_int(B_L, 4)));
        if(aux_l) {
          C.reserve(4, 1);
          C.tree.push_back(1);
          C.tree.push_back(1);
          C.tree.push_back(0);
          C.tree.push_back(0);
          C.l.push_back(aux_l);
        } else {
          C.reserve(2, 0);
          C.tree.push_back(1);
          C.tree.push_back(0);
          C.tree[0] = 1;
        }

        C.height_tree = curr_h;
        C.m = m;
        C.msize = msize;
        C.rmsize = rmsize;

        A_tree += 3;
        B_tree += 3;
        A_L += 4;
        B_L += 4;
        return;
      }

      uint64_t A_save = A_tree;
      // see if there is a pointer A
      if(P.size() > 0 && A_tree < tree.size() - 4 &&
         tree[A_tree + 1] && !tree[A_tree + 2] && !tree[A_tree + 3]) {
#ifdef DEBUG
        cout << "POINTER IN A" << endl;
#endif
        //uint64_t read_PoL = (A_tree >= occ_PoL.size() ? PoL.size() : rank1_occ_PoL(A_tree));
        uint64_t read_PoL = rank_occ_PoL(A_tree);
        if(!(read_PoL >= PoL.size() || access_PoL(read_PoL) == 0)) {
          uint64_t read_P = rank1_PoL(read_PoL);

          uint64_t where_to_move = select_real_tree(P[read_P] + 1);

          //uint64_t rank_p_here = rank1_occ_PoL(A_tree);
          uint64_t rank_p_here = rank_occ_PoL(A_tree);
          //uint64_t rank_p_to_move = rank1_occ_PoL(where_to_move);
          uint64_t rank_p_to_move = rank_occ_PoL(where_to_move);

          uint64_t pointer_A = rank1_PoL(rank_p_here);
          A_lvs_sk += (rank0_PoL(rank_p_here) - rank0_PoL(rank_p_to_move)) * 4;
          if(pointer_A > 0) {
            uint64_t pointer_A2 = rank1_PoL(rank_p_to_move);
            A_lvs_sk += pre_skips_A[pointer_A - 1] - (pointer_A2 > 0 ? pre_skips_A[pointer_A2 - 1] : 0);
          }
          A_tree = where_to_move;
        }
      }
      uint64_t B_save = B_tree;
      // see if there is a pointer B
      if(B.P.size() > 0 && B_tree < B.tree.size() - 4 &&
         B.tree[B_tree + 1] && !B.tree[B_tree + 2] && !B.tree[B_tree + 3]) {
#ifdef DEBUG
        cout << "POINTER IN B" << endl;
#endif
        //uint64_t read_PoL = (B_tree >= B.occ_PoL.size() ? B.PoL.size() : B.rank1_occ_PoL(B_tree));
        uint64_t read_PoL = B.rank_occ_PoL(B_tree);
        if(!(read_PoL >= B.PoL.size() || B.PoL[read_PoL] == 0)) {
          uint64_t read_P = B.rank1_PoL(read_PoL);

          uint64_t where_to_move = B.select_real_tree(B.P[read_P] + 1);
          //uint64_t rank_p_here = B.rank1_occ_PoL(B_tree);
          uint64_t rank_p_here = B.rank_occ_PoL(B_tree);
          //uint64_t rank_p_to_move = B.rank1_occ_PoL(where_to_move);
          uint64_t rank_p_to_move = B.rank_occ_PoL(where_to_move);

          uint64_t pointer_B = B.rank1_PoL(rank_p_here);
          B_lvs_sk += (B.rank0_PoL(rank_p_here) - B.rank0_PoL(rank_p_to_move)) * 4;
          if(pointer_B > 0) {
            uint64_t pointer_B2 = B.rank1_PoL(rank_p_to_move);
            B_lvs_sk += pre_skips_B[pointer_B - 1] - (pointer_B2 ? pre_skips_B[pointer_B2 - 1] : 0);
          }
          B_tree = where_to_move;
        }
      }
     // check if you need to comeback A
      //  A_0 | A_1
      //  ---------
      //  A_2 | A_3
      uint64_t A_0, A_1, A_2, A_3;
      uint64_t A_0_L, A_1_L, A_2_L, A_3_L;

      //  B_0 | B_1
      //  ---------
      //  B_2 | B_3
      uint64_t B_0, B_1, B_2, B_3;
      uint64_t B_0_L, B_1_L, B_2_L, B_3_L;

      //  C_0 | C_1
      //  ---------
      //  C_2 | C_3
      plain_tree C_0, C_1, C_2, C_3;
      plain_tree C_0_0, C_1_2, C_0_1, C_1_3, C_2_0, C_3_2, C_2_1, C_3_3;

      A_tree++;
      A_0 = A_tree;
      A_0_L = A_L;

      B_tree++;
      B_0 = B_tree;
      B_0_L = B_L;
      // A_0 * B_0
      mul(A_tree, A_L, A_L_S, pre_skips_A, A_lvs_sk, 0,
          B, B_tree, B_L, B_L_S, pre_skips_B, B_lvs_sk, 0,
          C_0_0, curr_h - 1);

      A_tree++;
      A_1 = A_tree;
      A_1_L = A_L;

      B_tree++;
      B_1 = B_tree;
      B_1_L = B_L;
      // A_0 * B_1
      mul(A_0, A_0_L, A_L_S, pre_skips_A, A_lvs_sk, 1,
          B, B_tree, B_L, B_L_S, pre_skips_B, B_lvs_sk, 0,
          C_0_1, curr_h - 1);

      B_tree++;
      B_2 = B_tree;
      B_2_L = B_L;
      // A_1 * B_2
      mul(A_tree, A_L, A_L_S, pre_skips_A, A_lvs_sk, 0,
          B, B_tree, B_L, B_L_S, pre_skips_B, B_lvs_sk, 0,
          C_1_2, curr_h - 1);

      C_0.reserve(2 * max(C_0_0.tree.size(), C_1_2.tree.size()), 2 * max(C_0_0.l.size(), C_1_2.l.size()));
      C_0_0.binsum(C_1_2, C_0);
      C_0_0.destroy();
      C_1_2.destroy();

      A_tree++;
      A_2 = A_tree;
      A_2_L = A_L;

      B_tree++;
      B_3 = B_tree;
      B_3_L = B_L;

      // A_1 * B_3
      mul(A_1, A_1_L, A_L_S, pre_skips_A, A_lvs_sk, 1,
          B, B_tree, B_L, B_L_S, pre_skips_B, B_lvs_sk, 0,
          C_1_3, curr_h - 1);

      C_1.reserve(2 * max(C_0_1.tree.size(), C_1_3.tree.size()), 2 * max(C_0_1.l.size(), C_1_3.l.size()));
      C_0_1.binsum(C_1_3, C_1);
      C_0_1.destroy();
      C_1_3.destroy();

      // A_2 * B_0
      mul(A_tree, A_L, A_L_S, pre_skips_A, A_lvs_sk, 0,
          B, B_0, B_0_L, B_L_S, pre_skips_B, B_lvs_sk, 1,
          C_2_0, curr_h - 1);

      A_tree++;
      A_3 = A_tree;
      A_3_L = A_L;

      // A_2 * B_1
      mul(A_2, A_2_L, A_L_S, pre_skips_A, A_lvs_sk, 1,
          B, B_1, B_1_L, B_L_S, pre_skips_B, B_lvs_sk, 1,
          C_2_1, curr_h - 1);

      // A_3 * B_2
      mul(A_tree, A_L, A_L_S, pre_skips_A, A_lvs_sk, 0,
          B, B_2, B_2_L, B_L_S, pre_skips_B, B_lvs_sk, 1,
          C_3_2, curr_h - 1);

      C_2.reserve(2 * max(C_2_0.tree.size(), C_3_2.tree.size()), 2 * max(C_2_0.l.size(), C_3_2.l.size()));
      C_2_0.binsum(C_3_2, C_2);
      C_2_0.destroy();
      C_3_2.destroy();

      // A_3 * B_3
      mul(A_3, A_3_L, A_L_S, pre_skips_A, A_lvs_sk, 1,
          B, B_3, B_3_L, B_L_S, pre_skips_B, B_lvs_sk, 1,
          C_3_3, curr_h - 1);

      C_3.reserve(2 * max(C_2_1.tree.size(), C_3_3.tree.size()), 2 * max(C_2_1.l.size(), C_3_3.l.size()));
      C_2_1.binsum(C_3_3, C_3);
      C_2_1.destroy();
      C_3_3.destroy();

      // merge results
#ifdef DEBUG
      cout << "MERGE" << endl;
      cout << C_0 << endl;
      cout << C_1 << endl;
      cout << C_2 << endl;
      cout << C_3 << endl;
#endif
      A_tree++;
      B_tree++;

      if(A_save > A_tree) A_tree = A_save + 3;
      if(B_save > B_tree) B_tree = B_save + 3;
      if(C_0.tree.size() == 2 &&
         C_1.tree.size() == 2 &&
         C_2.tree.size() == 2 &&
         C_3.tree.size() == 2) {
        C.tree.push_back(1);
        C.tree.push_back(0);
        C.height_tree = curr_h;
        C.m = m;
        C.msize = msize;
        C.rmsize = rmsize;
        return;
      }
      C.tree.reserve(2 + C_0.tree.size() + C_1.tree.size() + C_2.tree.size() + C_3.tree.size());
      C.tree.push_back(1);
      C.tree.insert(C.tree.end(), C_0.tree.begin(), C_0.tree.end());
      C.tree.insert(C.tree.end(), C_1.tree.begin(), C_1.tree.end());
      C.tree.insert(C.tree.end(), C_2.tree.begin(), C_2.tree.end());
      C.tree.insert(C.tree.end(), C_3.tree.begin(), C_3.tree.end());
      C.tree.push_back(0);

      C.l.reserve(C_0.l.size() + C_1.l.size() + C_2.l.size() + C_3.l.size());
      C.l.insert(C.l.end(), C_0.l.begin(), C_0.l.end());
      C.l.insert(C.l.end(), C_1.l.begin(), C_1.l.end());
      C.l.insert(C.l.end(), C_2.l.begin(), C_2.l.end());
      C.l.insert(C.l.end(), C_3.l.begin(), C_3.l.end());

      C.height_tree = curr_h;
      C.m = m;
      C.msize = msize;
      C.rmsize = rmsize;


      return;
    }

    void write(ofstream& out) {
      // writing integers first
      out.write((char*) &msize, sizeof(uint64_t));
      out.write((char*) &rmsize, sizeof(uint64_t));
      out.write((char*) &m, sizeof(uint64_t));
      out.write((char*) &height_tree, sizeof(uint64_t));
      out.write((char*) &last_bit_t, sizeof(uint64_t));
      out.write((char*) &last_bit_l, sizeof(uint64_t));
      P.serialize(out);

      //occ_PoL.serialize(out);
      //rank1_occ_PoL.serialize(out);

      PoL.serialize(out);

      rank1_PoL.serialize(out);
      rank0_PoL.serialize(out);
      select1_PoL.serialize(out);
      select0_PoL.serialize(out);

      real_tree.serialize(out);
      select_real_tree.serialize(out);

      tree.serialize(out);
      tree_support.serialize(out);
      l.serialize(out);
      blocks.serialize(out);
    }

    void load(ifstream& in) {
      // writing integers first
      in.read((char*) &msize, sizeof(uint64_t));
      in.read((char*) &rmsize, sizeof(uint64_t));
      in.read((char*) &m, sizeof(uint64_t));
      in.read((char*) &height_tree, sizeof(uint64_t));
      in.read((char*) &last_bit_t, sizeof(uint64_t));
      in.read((char*) &last_bit_l, sizeof(uint64_t));

      P.load(in);

      //occ_PoL.load(in);
      //rank1_occ_PoL.load(in, &occ_PoL);

      PoL.load(in);
      rank1_PoL.load(in, &PoL);
      rank0_PoL.load(in, &PoL);
      select1_PoL.load(in, &PoL);
      select0_PoL.load(in, &PoL);

      real_tree.load(in);
      select_real_tree.load(in, &real_tree);

      tree.load(in);
      tree_support.load(in, &tree);

      l.load(in);
      blocks.load(in);
    }

    uint64_t size_in_bits() {
      uint64_t total = sizeof(uint64_t) * 4 +
             size_in_bytes(tree) * 8 +
             size_in_bytes(tree_support) * 8 +
             size_in_bytes(l) * 8 +
             size_in_bytes(P) * 8 +
             (size_in_bytes(real_tree) + size_in_bytes(select_real_tree)) * 8 +
             size_in_bytes(PoL) * 8 + size_in_bytes(rank1_PoL) * 8 + size_in_bytes(rank0_PoL) * 8 +
             size_in_bytes(select1_PoL) * 8 + size_in_bytes(select0_PoL) * 8 + 
             size_in_bytes(blocks) * 8;
             //size_in_bytes(occ_PoL) * 8 + size_in_bytes(rank1_occ_PoL) * 8;
#ifdef INFO_SPACE
      cout << "BITS" << endl;
      cout << "  Tree        : " << (size_in_bytes(tree)) * 8 << " " << (double) 100 *  (size_in_bytes(tree)) * 8 / total << endl;
      cout << "  Tree Support: " << (size_in_bytes(tree_support)) * 8 << " " << (double) 100 *  (size_in_bytes(tree_support)) * 8 / total  << endl;
      cout << "  L           : " << (size_in_bytes(l)) * 8 << " " << (double) 100 *  (size_in_bytes(l)) * 8 / total << endl;
      cout << "  P           : " << (size_in_bytes(P)) * 8 << " " << (double) 100 *  (size_in_bytes(P)) * 8 / total << endl;
      cout << "  Real P      : " << (size_in_bytes(real_tree) + size_in_bytes(select_real_tree)) * 8 << " " << (double) 100 *  (size_in_bytes(real_tree) + size_in_bytes(select_real_tree)) * 8 / total << endl;
      cout << "  PoL         : " << (size_in_bytes(PoL) + size_in_bytes(rank1_PoL) + size_in_bytes(rank0_PoL) + size_in_bytes(select1_PoL) + size_in_bytes(select0_PoL)) * 8 << " " << (double) 100 *  (size_in_bytes(PoL) + size_in_bytes(rank1_PoL) + size_in_bytes(rank0_PoL) + size_in_bytes(select1_PoL) + size_in_bytes(select0_PoL)) * 8 / total << endl;
      //cout << "  occ_PoL     : " << (size_in_bytes(occ_PoL) + size_in_bytes(rank1_occ_PoL)) * 8 << " " << (double) 100 *  (size_in_bytes(occ_PoL) + size_in_bytes(rank1_occ_PoL)) * 8 / total << endl;
      cout << "  occ_PoL (bs): " << (size_in_bytes(blocks)) * 8 << " " << (double) 100 * (size_in_bytes(blocks)) * 8 / total << endl;
#endif
      return total;
    }

    friend ostream& operator<<(ostream& os, const k2tree_bp_sdsl_idems< k, bv_leaves, bit_vector_1, rank1_1, bit_vector_2, rank1_2, rank0_2, select1_2, select0_2, b_size > &k2tree) {
      cout << "HT  : " << k2tree.height_tree << endl;
      cout << "Tree: ";
      for(uint64_t i = 0; i < k2tree.tree.size(); i++) {
        cout << (k2tree.tree[i] ? "(" : ")");
      }
      cout << endl;
      cout << "L   : ";
      for(uint64_t i = 0; i < k2tree.l.size(); i++) {
        if(i % 4 == 0 && !(i == 0)) cout << " ";
        cout << (k2tree.l[i] ? "1" : "0");
      }
      cout << endl;
      cout << "P   : ";
      for(uint64_t i = 0; i < k2tree.P.size(); i++) {
        cout << k2tree.P[i] << " ";
      }
      cout << endl;
      cout << "RT  : ";
      for(uint64_t i = 0; i < k2tree.real_tree.size(); i++) {
        cout << k2tree.real_tree[i];
      }
      cout << endl;
      cout << "PoL : ";
      for(uint64_t i = 0; i < k2tree.PoL.size(); i++) {
        cout << k2tree.PoL[i];
      }
      cout << endl;
//      cout << "OPoL: ";
//      for(uint64_t i = 0; i < k2tree.occ_PoL.size(); i++) {
//        cout << k2tree.occ_PoL[i];
//      }
//      cout << endl;
      cout << "OPoL: ";
      for(uint64_t i = 0; i < k2tree.blocks.size(); i++) {
        cout << k2tree.blocks[i] << " ";
      }
      return os;
    }
};
#endif //K2_TREE_BP_SDSL_IDEM
