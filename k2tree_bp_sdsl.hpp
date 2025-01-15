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
#include <sdsl/rrr_vector.hpp>
#include <sdsl/sd_vector.hpp>
#include <stack>
#include <type_traits>
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

    rrr_vector<127> l; // real values
    uint64_t last_bit_l; // universe

    sd_vector<> leaves;
    rank_support_sd<> rank_leaves;

    uint64_t msize;
    uint64_t rmsize;
    uint64_t m;

    void add_one(vector< uint64_t > &bv, uint64_t &pos_to_add) {
      bv.push_back(pos_to_add++);
    }

    void add_zero(vector< uint64_t > &bv, uint64_t &pos_to_add) {
      pos_to_add++;
    }

    uint64_t bits_interleave(int64_t a, int64_t b) {
      uint64_t r = 0;
      assert(a<=UINT32_MAX && b <= UINT32_MAX);
      int c = 0;
      while(a!=0 || b!=0) {
        r |= (b&1)<<c++;
        r |= (a&1)<<c++;
        a >>= 1; b>>=1;  
        assert(c<=64);
      }
      return r;
    }

    size_t binsearch(vector< uint64_t >::iterator ia, size_t n, uint64_t x) {
      assert(n>0);
      size_t l=0, r=n-1;
      while(l<r) {
        size_t m = (l+r)/2;
        if(ia[m]<x) l=m+1;
        else if(ia[m]==x) return m;
        else r=m; // replace with r = m-1 and later return r+1?
      }
      assert(l==r);
      if(ia[l]<x) {
        assert(r==n-1);
        return n;   // replace with return r+1?
      }
      return l;
    }

  public:
    uint64_t nodes() { return (tree_support.find_close(0) + 1) / 2; }

    k2tree_bp_sdsl() {}

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
      
      vector< uint64_t > ia_ones;
      for(const auto& one : ones) {
        ia_ones.push_back(bits_interleave(one.first, one.second));
      }

      sort(ia_ones.begin(), ia_ones.end());

      // first i will do it asuming k = 2
      // then i will generalize
#ifdef DEBUG
      string balance_string = "";
      cout << "Initialize recursion..." << endl;
#endif // DEBUG
      stack< tuple< uint64_t, uint64_t, uint64_t, uint64_t, vector< uint64_t >::iterator, uint64_t, bool, bool > > recursion;
      recursion.push(make_tuple(msize, 0, 0, 0, ia_ones.begin(), ia_ones.size(), true, false));

      vector< uint64_t > bv_tree;
      vector< uint64_t > bv_l;

      uint64_t pos_to_add = 0;
      uint64_t pos_to_add_l = 0;

      while(!recursion.empty()) {
        auto [subm_size, init_x, init_y, smin, ia, n_ia, one_one, flag] = recursion.top();
#ifdef DEBUG
        cout << "Recursion call..." << endl;
        cout << "Current Sub Matrix Size: " << subm_size << " init x: " << init_x << " init y: " << init_y << " visited: " << flag << endl;
        cout << " n ia: " << n_ia << endl;
#endif // DEBUG
        recursion.pop();
        
        if(flag) {
#ifdef DEBUG
          balance_string += ")";
          cout << "Adding )..." << endl;
#endif // DEBUG
          add_zero(bv_tree, pos_to_add);
          continue;
        }

        recursion.push(make_tuple(subm_size, init_x, init_y, smin, ia, n_ia, one_one, true));
#ifdef DEBUG
        balance_string += "(";
        cout << "Adding (..." << endl;
#endif // DEBUG
        add_one(bv_tree, pos_to_add);

        if(!one_one) {
          continue;
        }

        
        if(subm_size == k) {
          vector< int64_t > t(4, 0);

          for(size_t i = 0; i < n_ia; i++) {
#ifdef DEBUG
            cout << "IA[" << i << "] = " << ia[i] << endl;
            cout << "smin = " << smin << endl;
#endif
            int64_t pos = (int64_t) (ia[i] - smin);
#ifdef DEBUG
            cout << "pos = " << pos << endl;
#endif
            t[pos] = 1;
          }

          for(const auto& bit : t) {
            if(bit) add_one(bv_l, pos_to_add_l);
            else add_zero(bv_l, pos_to_add_l);
          }

#ifdef DEBUG
          balance_string += "(";
          cout << "Adding (..." << endl;
#endif // DEBUG
          add_one(bv_tree, pos_to_add);
#ifdef DEBUG
          balance_string += ")";
          cout << "Adding )..." << endl;
#endif // DEBUG
          add_zero(bv_tree, pos_to_add);
          continue;
        }

        uint64_t range = (subm_size / 2) * (subm_size / 2);
        uint64_t left = smin + range;
        uint64_t mid = left + range;
        uint64_t right = mid + range;

        size_t imid = binsearch(ia, n_ia, mid);
        size_t ileft = imid >0 ? binsearch(ia, imid, left) : 0;
        size_t iright = imid < n_ia ? binsearch(ia + imid, n_ia - imid, right) + imid : n_ia;

        if(iright < n_ia) { // right-bot 
          recursion.push(make_tuple(subm_size / 2, init_x + subm_size / 2, init_y + subm_size / 2, right, ia + iright, n_ia - iright, true, false));
        } else {
          recursion.push(make_tuple(subm_size / 2, init_x + subm_size / 2, init_y + subm_size / 2, right, ia + iright, n_ia - iright, false, false));
        }
        
        if(iright > imid) { // left-bot 
          recursion.push(make_tuple(subm_size / 2, init_x + subm_size / 2, init_y, mid, ia + imid, iright - imid, true, false));
        } else {
          recursion.push(make_tuple(subm_size / 2, init_x + subm_size / 2, init_y, mid, ia + imid, iright - imid, false, false));
        }
        
        if(ileft < imid) { // right-up
          recursion.push(make_tuple(subm_size / 2, init_x, init_y + subm_size / 2, left, ia + ileft, imid - ileft, true, false));
        } else {
          recursion.push(make_tuple(subm_size / 2, init_x, init_y + subm_size / 2, left, ia + ileft, imid - ileft, false, false));
        }

        if(ileft > 0) { // left-up
          recursion.push(make_tuple(subm_size / 2, init_x, init_y, smin, ia, ileft, true, false));
        } else {
          recursion.push(make_tuple(subm_size / 2, init_x, init_y, smin, ia, ileft, false, false));
        }
      }

      last_bit_t = pos_to_add;
      last_bit_l = pos_to_add_l;

#ifdef DEBUG
      cout << "Result: " << balance_string << "..." << endl;
      cout << "Init tree " << pos_to_add << "..." << endl;
#endif // DEBUG
      tree = bit_vector(pos_to_add, 0);
      for(const auto& bit : bv_tree) tree[bit] = 1;

      // i think this can be improved
      bit_vector aux_leaves(pos_to_add, 0);
      for(uint64_t i = 0; i < pos_to_add - 3; i++) {
        if(tree[i] && tree[i + 1] && !tree[i + 2] && !tree[i + 3])
          aux_leaves[i] = 1;
      }

      leaves = sd_vector<>(aux_leaves);
      util::init_support(rank_leaves, &leaves);

#ifdef DEBUG
      cout << "Init L " << pos_to_add_l << "..." << endl;
#endif // DEBUG
      auto aux_l = bit_vector(pos_to_add_l, 0);
      for(const auto& bit : bv_l) aux_l[bit] = 1;
      l = rrr_vector<127>(aux_l);

#ifdef DEBUG
      cout << "Init Tree support..." << endl;
#endif // DEBUG
      tree_support = bp_support_sada<>(&tree);

#ifdef DEBUG
      cout << "End k2tree building..." << endl;
#endif // DEBUG
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

    k2tree_bp_sdsl<k> operator|(const k2tree_bp_sdsl<k>& B) {
      uint64_t pa, pb;
      uint64_t pLa, pLb;

      pa = pb = pLa = pLb = 0;

      k2tree_bp_sdsl<k> C;

      if(B.tree.size() == 2 && tree.size() == 2) {
        C.tree = bit_vector(2, 0);
        C.tree[0] = 1;
        C.height_tree = height_tree;
        C.msize = msize;
        C.rmsize = rmsize;
        C.m = m;
      }

      if(B.tree.size() == 2) {
        C.tree = tree;
        C.l = l;
        C.tree_support = tree_support;
        C.height_tree = height_tree;
        C.msize = msize;
        C.rmsize = rmsize;
        C.m = m;
        return C;
      }

      if(tree.size() == 2) {
        C.tree = B.tree;
        C.l = B.l;
        C.tree_support = B.tree_support;
        C.height_tree = B.height_tree;
        C.msize = B.msize;
        C.rmsize = B.rmsize;
        C.m = B.m;
        return C;
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

      auto aux_l = bit_vector(curr_bit_L, 0);
      for(auto bit : bits_L) aux_l[bit] = 1;
      C.l = rrr_vector<127>(aux_l);

      C.tree_support = bp_support_sada<>(&C.tree);

      C.last_bit_t = curr_bit_tree;
      C.last_bit_l = curr_bit_L;
      C.height_tree = height_tree;
      C.msize = msize;
      C.rmsize = rmsize;
      C.m = m;
 
      // i think this can be improved
      bit_vector aux_leaves(curr_bit_tree, 0);
      for(uint64_t i = 0; i < curr_bit_tree - 3; i++) {
        if(C.tree[i] && C.tree[i + 1] && !C.tree[i + 2] && !C.tree[i + 3])
          aux_leaves[i] = 1;
      }

      C.leaves = sd_vector<>(aux_leaves);
      util::init_support(C.rank_leaves, &C.leaves);     
      return C;
    }

    k2tree_bp_sdsl<k> operator*(const k2tree_bp_sdsl<k>& B) {
      return mul(0, 0, B, 0, 0, height_tree);
    }   

    k2tree_bp_sdsl<k> mul(uint64_t A_tree, uint64_t A_L,
                          const k2tree_bp_sdsl<k> &B, uint64_t B_tree, uint64_t B_L,
                          uint64_t curr_h) {
#ifdef DEBUG
      cout << "Current Height: " << curr_h << endl;
      cout << "Start Matrix A: " << A_tree << endl;
      cout << "Start Matrix B: " << B_tree << endl;
#endif
      k2tree_bp_sdsl<k> C;
      // submatrix A or B full of 0's
      if((tree[A_tree] && !tree[A_tree + 1]) ||
         (B.tree[B_tree] && !B.tree[B_tree + 1])) { 
#ifdef DEBUG
        cout << "Full of 0's" << endl;
#endif
        C.tree = bit_vector(2, 0);
        C.tree[0] = 1;
        C.height_tree = curr_h;
        C.m = m;
        C.msize = msize;
        C.rmsize = rmsize;
        return C;
      }

      // base case, leave 
      if(curr_h == 1) { 
#ifdef DEBUG
        cout << "Leaf!" << endl;
#endif
        auto aux_l = bit_vector(4, 0);
        bool found_1 = 0;
        aux_l[0] = (l[A_L] & B.l[B_L]) | (l[A_L + 1] * B.l[B_L + 2]);
        found_1 = found_1 | aux_l[0];
        aux_l[1] = (l[A_L] & B.l[B_L + 1]) | (l[A_L + 1] * B.l[B_L + 3]);
        found_1 = found_1 | aux_l[1];
        aux_l[2] = (l[A_L + 2] & B.l[B_L]) | (l[A_L + 3] * B.l[B_L + 2]);
        found_1 = found_1 | aux_l[2];
        aux_l[3] = (l[A_L + 2] & B.l[B_L + 1]) | (l[A_L + 3] * B.l[B_L + 3]);
        found_1 = found_1 | aux_l[3];

        if(found_1) {
          C.tree = bit_vector(4, 0);
          C.tree[0] = 1;
          C.tree[1] = 1;
          C.l = rrr_vector<127>(aux_l);
        } else {
          C.tree = bit_vector(2, 0);
          C.tree[0] = 1;
        }
        C.height_tree = curr_h;
        C.m = m;
        C.msize = msize;
        C.rmsize = rmsize;
        return C;
      }

 

      //  A_0 | A_1
      //  ---------
      //  A_2 | A_3
      uint64_t A_0 = A_tree + 1;
      uint64_t A_L_0 = rank_leaves(A_0) * 4;
#ifdef DEBUG
      cout << "A_0" << endl;
      cout << "  " << A_0 << endl;
      cout << "  " << A_L_0 << endl;
#endif
      uint64_t A_1 = tree_support.find_close(A_0) + 1;
      uint64_t A_L_1 = rank_leaves(A_1) * 4;
#ifdef DEBUG
      cout << "A_1" << endl;
      cout << "  " << A_1 << endl;
      cout << "  " << A_L_1 << endl;
#endif
      uint64_t A_2 = tree_support.find_close(A_1) + 1;
      uint64_t A_L_2 = rank_leaves(A_2) * 4;
#ifdef DEBUG
      cout << "A_2" << endl;
      cout << "  " << A_2 << endl;
      cout << "  " << A_L_2 << endl;
#endif
      uint64_t A_3 = tree_support.find_close(A_2) + 1;
      uint64_t A_L_3 = rank_leaves(A_3) * 4;
#ifdef DEBUG
      cout << "A_3" << endl;
      cout << "  " << A_3 << endl;
      cout << "  " << A_L_3 << endl;
#endif

      //  B_0 | B_1
      //  ---------
      //  B_2 | B_3
      uint64_t B_0 = B_tree + 1;
      uint64_t B_L_0 = B.rank_leaves(B_0) * 4;
#ifdef DEBUG
      cout << "B_0" << endl;
      cout << "  " << B_0 << endl;
      cout << "  " << B_L_0 << endl;
#endif
      uint64_t B_1 = B.tree_support.find_close(B_0) + 1;
      uint64_t B_L_1 = B.rank_leaves(B_1) * 4;
#ifdef DEBUG
      cout << "B_1" << endl;
      cout << "  " << B_1 << endl;
      cout << "  " << B_L_1 << endl;
#endif
      uint64_t B_2 = B.tree_support.find_close(B_1) + 1;
      uint64_t B_L_2 = B.rank_leaves(B_2) * 4;
#ifdef DEBUG
      cout << "B_2" << endl;
      cout << "  " << B_2 << endl;
      cout << "  " << B_L_2 << endl;
#endif
      uint64_t B_3 = B.tree_support.find_close(B_2) + 1;
      uint64_t B_L_3 = B.rank_leaves(B_3) * 4;
#ifdef DEBUG
      cout << "B_3" << endl;
      cout << "  " << B_3 << endl;
      cout << "  " << B_L_3 << endl;
#endif

      //  C_0 | C_1
      //  ---------
      //  C_2 | C_3
      auto C_0 = mul(A_0, A_L_0, B, B_0, B_L_0, curr_h - 1) | mul(A_1, A_L_1, B, B_2, B_L_2, curr_h - 1);
      auto C_1 = mul(A_0, A_L_0, B, B_1, B_L_1, curr_h - 1) | mul(A_1, A_L_1, B, B_3, B_L_3, curr_h - 1);
      auto C_2 = mul(A_2, A_L_2, B, B_0, B_L_0, curr_h - 1) | mul(A_3, A_L_3, B, B_2, B_L_2, curr_h - 1);
      auto C_3 = mul(A_2, A_L_2, B, B_1, B_L_1, curr_h - 1) | mul(A_3, A_L_3, B, B_3, B_L_3, curr_h - 1);

      if(C_0.tree.size() == 2 &&
         C_1.tree.size() == 2 &&
         C_2.tree.size() == 2 &&
         C_3.tree.size() == 2) {
        C.tree = bit_vector(2, 0);
        C.tree[0] = 1;
        C.height_tree = curr_h;
        C.m = m;
        C.msize = msize;
        C.rmsize = rmsize;
        return C;
      }

      // merge results
#ifdef DEBUG
      cout << "MERGE" << endl;
      cout << C_0 << endl;
      cout << C_1 << endl;
      cout << C_2 << endl;
      cout << C_3 << endl;
#endif
      C.tree = bit_vector(C_0.tree.size() + C_1.tree.size() + C_2.tree.size() + C_3.tree.size() + 2, 0);
      C.tree[0] = 1;
      for(uint64_t i = 0; i < C_0.tree.size(); i++)
        C.tree[1 + i] = C_0.tree[i];
      for(uint64_t i = 0; i < C_1.tree.size(); i++)
        C.tree[1 + C_0.tree.size() + i] = C_1.tree[i];
      for(uint64_t i = 0; i < C_2.tree.size(); i++)
        C.tree[1 + C_0.tree.size() + C_1.tree.size() + i] = C_2.tree[i];
      for(uint64_t i = 0; i < C_3.tree.size(); i++)
        C.tree[1 + C_0.tree.size() + C_1.tree.size() + C_2.tree.size() + i] = C_3.tree[i];

      auto aux_C_L = bit_vector(C_0.l.size() + C_1.l.size() + C_2.l.size() + C_3.l.size(), 0);
      for(uint64_t i = 0; i < C_0.l.size(); i++)
        aux_C_L[i] = C_0.l[i];
      for(uint64_t i = 0; i < C_1.l.size(); i++)
        aux_C_L[C_0.l.size() + i] = C_1.l[i];
      for(uint64_t i = 0; i < C_2.l.size(); i++)
        aux_C_L[C_0.l.size() + C_1.l.size() + i] = C_2.l[i];
      for(uint64_t i = 0; i < C_3.l.size(); i++)
        aux_C_L[C_0.l.size() + C_1.l.size() + C_2.l.size() + i] = C_3.l[i];

      C.l = rrr_vector<127>(aux_C_L);
      C.height_tree = curr_h;
      C.m = m;
      C.msize = msize;
      C.rmsize = rmsize;


      return C;
    }

    uint64_t size_in_bits() {
      return sizeof(uint64_t) * 5 +
             size_in_bytes(tree) * 8 +
             size_in_bytes(tree_support) * 8 +
             size_in_bytes(l) * 8 +
             size_in_bytes(leaves) * 8 + size_in_bytes(rank_leaves) * 8;
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
      cout << endl;
      cout << "Leaves: ";
      for(uint64_t i = 0; i < k2tree.leaves.size(); i++) {
        cout << (k2tree.leaves[i] ? "1" : "0");
      }
      return os;
    }
};
#endif // !K2_TREE_BP_SDSL
