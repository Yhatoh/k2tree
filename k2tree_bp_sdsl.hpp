#ifndef K2_TREE_BP_SDSL
#define K2_TREE_BP_SDSL

// std includes
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <stack>
#include <utility>
#include <vector>
#include <tuple>

// sdsl includes
#include <sdsl/int_vector.hpp>
#include <sdsl/bp_support_sada.hpp>
#include <sdsl/util.hpp>
#include <sdsl/lcp.hpp>
#include <sdsl/construct.hpp>
#include <sdsl/io.hpp>
#include <sdsl/lcp_support_sada.hpp>
#include <sdsl/rrr_vector.hpp>
#include <sdsl/sd_vector.hpp>

// local includes
#include "util.hpp"
#include "plaintree.hpp"

using namespace std;
using namespace sdsl;

// k2-tree
// parameters:
//   * k * k: amount of children per node
template< uint64_t k = 2, class bv_leaves = bit_vector >
class k2tree_bp_sdsl {
  public:
    uint64_t height_tree;

    bp_support_sada<> tree_support;
    bit_vector tree; // k2tree
    uint64_t last_bit_t; // universe

    bv_leaves l; // real values
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
    uint64_t size() {
      uint64_t m = 0;
      sdsl::rank_support_rrr<1, 127> rank(&l);
      return rank(l.size());
    }

    uint64_t size_matrix() { return rmsize; }
    uint64_t nodes() { return tree.size() / 2; }

    k2tree_bp_sdsl() {}
    
    k2tree_bp_sdsl(plain_tree &pd) {
      tree = bit_vector(pd.tree.size(), 0);
      for(uint64_t i = 0; i < pd.tree.size(); i++) tree[i] = pd.tree[i];

      last_bit_t = tree.size();
      tree_support = bp_support_sada<>(&tree);

      //l = pd.l;
      bit_vector aux_l = bit_vector(pd.l.size() * 4, 0);
      for(uint64_t i = 0; i < pd.l.size(); i++) {
        for(uint64_t j = 0; j < 4; j++) {
          if(pd.l[i] & (1 << j)) aux_l[i * 4 + j] = 1;
        }
      }
      l = bv_leaves(aux_l);
      last_bit_l = l.size();

      height_tree = pd.height_tree;
      msize = pd.msize;
      rmsize = pd.rmsize;
      m = pd.m;

      // i think this can be improved
      bit_vector aux_leaves(tree.size(), 0);
      for(uint64_t i = 0; i < tree.size() - 3; i++) {
        if(tree[i] && tree[i + 1] && !tree[i + 2] && !tree[i + 3])
          aux_leaves[i] = 1;
      }

      leaves = sd_vector<>(aux_leaves);
      util::init_support(rank_leaves, &leaves);
    }


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
      l = bv_leaves(aux_l);

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

    /*
    void binsum(const k2tree_bp_sdsl<k, bv_leaves>& B, plain_tree &C) {
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
        uint8_t num = l[0];
        for(uint64_t A_L = 1; A_L < l.size(); A_L++) {
          if(A_L % 4 == 0) {
            C.l.push_back(num);
            num = l[A_L];
          } else {
            num |= l[A_L] << (A_L % 4);
          }
        }
        C.height_tree = height_tree;
        C.msize = msize;
        C.rmsize = rmsize;
        C.m = m;
        return;
      }

      if(tree.size() == 2) {
        C.tree = B.tree;
        uint8_t num = B.l[0];
        for(uint64_t B_L = 1; B_L < B.l.size(); B_L++) {
          if(B_L % 4 == 0) {
            C.l.push_back(num);
            num = B.l[B_L];
          } else {
            num |= B.l[B_L] << (B_L % 4);
          }
        }
        C.l.push_back(num);
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
            C.l |= (l[A_L] | B.l[B_L]) << 1;
            A_L++;
            B_L++;
          }
          A_tree++;
          B_tree++;
          curr_depth++;
        } else if(tree[A_tree] && !B.tree[B_tree]) {
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
                C.l |= (l[A_L]) << 1;
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
        } else if(!tree[A_tree] && B.tree[B_tree]) {
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
              for(uint64_t i = 0; i < 4; i++) {
                C.l |= (B.l[B_L]) << 1;
                B_L++;
              }
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

      C.height_tree = height_tree;
      C.msize = msize;
      C.rmsize = rmsize;
      C.m = m;
      return;
    }
*/

    void mul(const k2tree_bp_sdsl<k, bv_leaves> &B, plain_tree &C) {
      uint64_t A_tree, B_tree;
      A_tree = B_tree = 0;
      uint64_t A_L, B_L;
      A_L = B_L = 0;
      sdsl::int_vector<4> A_L_S(l.size() / 4, 0);
      sdsl::int_vector<4> B_L_S(B.l.size() / 4, 0);

      mul(A_tree, A_L, 0, A_L_S, B, B_tree, B_L, 0, B_L_S, C, height_tree);
    }

    void mul(uint64_t &A_tree, uint64_t &A_L, bool A_flag, sdsl::int_vector<4> &A_L_S,
             const k2tree_bp_sdsl<k, bv_leaves> &B, uint64_t &B_tree, uint64_t &B_L, bool B_flag, sdsl::int_vector<4> &B_L_S,
             plain_tree &C,
             uint8_t curr_h) {
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
        B_L = B.rank_leaves(B_tree) * 4;
        return;
      } else if(B_f0) {
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
        A_L = rank_leaves(A_tree) * 4;
        return;
      }

      // base case, leave 
      if(curr_h == 1) { 
#ifdef DEBUG
        cout << "Leaf!" << endl;
#endif
        uint8_t aux_l = minimat_mul((A_L_S[A_L >> 2] ? A_L_S[A_L >> 2] : A_L_S[A_L >> 2] = l.get_int(A_L, 4)),
                                    (B_L_S[B_L >> 2] ? B_L_S[B_L >> 2] : B_L_S[B_L >> 2] = B.l.get_int(B_L, 4)));

        if(aux_l > 0) {
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
      mul(A_tree, A_L, 0, A_L_S, B, B_tree, B_L, 0, B_L_S, C_0_0, curr_h - 1); // A_tree == A_1 && B_tree == B_1
      
      A_tree++;
      A_1 = A_tree;
      A_1_L = A_L;

      B_tree++;
      B_1 = B_tree;
      B_1_L = B_L;
      // A_0 * B_1
      mul(A_0, A_0_L, 1, A_L_S, B, B_tree, B_L, 0, B_L_S, C_0_1, curr_h - 1); // A_tree == A_1 && B_tree == B_2

      B_tree++;
      B_2 = B_tree;
      B_2_L = B_L;
      // A_1 * B_2
      mul(A_tree, A_L, 0, A_L_S, B, B_tree, B_L, 0, B_L_S, C_1_2, curr_h - 1);

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
      mul(A_1, A_1_L, 1, A_L_S, B, B_tree, B_L, 0, B_L_S, C_1_3, curr_h - 1);

      C_1.reserve(2 * max(C_0_1.tree.size(), C_1_3.tree.size()), 2 * max(C_0_1.l.size(), C_1_3.l.size()));
      C_0_1.binsum(C_1_3, C_1);
      C_0_1.destroy();
      C_1_3.destroy();

      // A_2 * B_0
      mul(A_tree, A_L, 0, A_L_S, B, B_0, B_0_L, 1, B_L_S, C_2_0, curr_h - 1);

      A_tree++;
      A_3 = A_tree;
      A_3_L = A_L;

      // A_2 * B_1
      mul(A_2, A_2_L, 1, A_L_S, B, B_1, B_1_L, 1, B_L_S, C_2_1, curr_h - 1);

      // A_3 * B_2
      mul(A_tree, A_L, 0, A_L_S, B, B_2, B_2_L, 1, B_L_S, C_3_2, curr_h - 1);

      C_2.reserve(2 * max(C_2_0.tree.size(), C_3_2.tree.size()), 2 * max(C_2_0.l.size(), C_3_2.l.size()));
      C_2_0.binsum(C_3_2, C_2);
      C_2_0.destroy();
      C_3_2.destroy();

      // A_3 * B_3
      mul(A_3, A_3_L, 1, A_L_S, B, B_3, B_3_L, 1, B_L_S, C_3_3, curr_h - 1);

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

      leaves.serialize(out);
      rank_leaves.serialize(out);

      tree.serialize(out);
      tree_support.serialize(out);
      l.serialize(out);
    }

    void load(ifstream& in) {
      // writing integers first
      in.read((char*) &msize, sizeof(uint64_t));
      in.read((char*) &rmsize, sizeof(uint64_t));
      in.read((char*) &m, sizeof(uint64_t));
      in.read((char*) &height_tree, sizeof(uint64_t));
      in.read((char*) &last_bit_t, sizeof(uint64_t));
      in.read((char*) &last_bit_l, sizeof(uint64_t));

      leaves.load(in);
      rank_leaves.load(in, &leaves);

      tree.load(in);
      tree_support.load(in, &tree);

      l.load(in);
    }

    uint64_t size_in_bits() {
      uint64_t total = sizeof(uint64_t) * 5 +
             size_in_bytes(tree) * 8 +
             size_in_bytes(tree_support) * 8 +
             size_in_bytes(l) * 8 +
             size_in_bytes(leaves) * 8 + size_in_bytes(rank_leaves) * 8;
#ifdef INFO_SPACE
      cout << "Leaves:" << rank_leaves(leaves.size()) << endl;
      cout << "BITS" << endl;
      cout << "  Tree        : " << (size_in_bytes(tree)) * 8 << " " << (double) (size_in_bytes(tree)) * 8 / total<< endl;
      cout << "  Tree Support: " << (size_in_bytes(tree_support)) * 8 << " " << (double) (size_in_bytes(tree_support)) * 8 / total << endl;
      cout << "  L           : " << (size_in_bytes(l)) * 8 << " " << (double) (size_in_bytes(l)) * 8 / total << endl;
      cout << "  leaves      : " << (size_in_bytes(leaves) + size_in_bytes(rank_leaves)) * 8  << " " << (double) (size_in_bytes(leaves) + size_in_bytes(rank_leaves)) * 8 / total<< endl;
#endif
      return total;
    }

    friend ostream& operator<<(ostream& os, const k2tree_bp_sdsl<k, bv_leaves> &k2tree) {
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
      cout << "Lvs : ";
      for(uint64_t i = 0; i < k2tree.leaves.size(); i++) {
        cout << (k2tree.leaves[i] ? "1" : "0");
      }
      return os;
    }
};
#endif // !K2_TREE_BP_SDSL
