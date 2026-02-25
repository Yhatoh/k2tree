#ifndef K2_TREE_BP_SDSL
#define K2_TREE_BP_SDSL

// std includes
#include <algorithm>
#include <cassert>
#include <cmath>
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
#include <sdsl/construct.hpp>
#include <sdsl/io.hpp>
#include <sdsl/rrr_vector.hpp>
#include <sdsl/sd_vector.hpp>

// local includes
#include "util.hpp"
#include "plaintree.hpp"
#include "debug.hpp"

using namespace std;
using namespace sdsl;

#define NUM_SUPPORT 36
#define GET_NODES(x) (x & ((1LL << NUM_SUPPORT) - 1))
#define GET_SKIPS(x) (x >> NUM_SUPPORT)
#define ENCODE(x, y) (x << NUM_SUPPORT) | y

struct child_info {
  uint64_t size_tree;
  uint32_t n_leaves;

  child_info(uint64_t x, uint64_t y, uint32_t leaves) {
    n_leaves = leaves;
    size_tree = ENCODE(x, y);
  }
  child_info() {}
};

struct traverse_info {
  uint64_t pos, l, node, size, n_l;
  bool info; // if have info or no
  
  traverse_info(uint64_t pos_, uint64_t l_,
                uint64_t node_, uint64_t size_,
                uint64_t n_l_, bool info_) :
    pos(pos_), l(l_), node(node_), size(size_), n_l(n_l_), info(info_) {}

  traverse_info() { pos = l = node = size = n_l = 0; }
};

// k2-tree
// parameters:
//   * k * k: amount of children per node
template< uint64_t k = 2, class bv_leaves = bit_vector >
class k2_bp {
  public:
    uint64_t height_tree;

    bp_support_sada<> tree_support;
    std::vector< child_info > child_support;
    std::vector< child_info > dynamic_support;
    bit_vector tree; // k2tree
    uint64_t last_bit_t; // universe

    bv_leaves l; // real values
    uint64_t last_bit_l; // universe

    sd_vector<> leaves;
    rank_support_sd<> rank_leaves;

    uint64_t msize;
    uint64_t rmsize;
    uint64_t m;

    uint64_t threshold;

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

    // x numbers to skip in child support
    // y amount of nodes

    void init_support_child(uint64_t &pos, std::vector< child_info > &c, uint64_t &leaves, uint64_t threshold) {
      assert(tree.size() > 0);
      assert(tree[pos] > 0); // should be always a (
      if(tree[pos] && tree[pos + 1] && !tree[pos + 2] && !tree[pos + 3]) {
        // jumping last level tree
        pos += 4;
        leaves++;
        return;
      }

      if(tree[pos] && !tree[pos + 1]) {
        // ignore leaf
        pos += 2;
        return ;
      }

      std::vector< std::vector< child_info > > each_child_size(4);
      std::vector< uint64_t > c_sizes(4, 0);
      std::vector< uint64_t > c_leaves(4, 0);

      pos++;
      uint64_t curr_pos = pos;
      init_support_child(pos, each_child_size[0], c_leaves[0], threshold);
      c_sizes[0] = (pos - curr_pos) / 2;

      //c.push_back(ENCODE(each_child_size[0].size(), (pos - curr_pos) / 2));

      curr_pos = pos;
      init_support_child(pos, each_child_size[1], c_leaves[1], threshold);
      c_sizes[1] = (pos - curr_pos) / 2;

      //c.push_back(ENCODE(each_child_size[1].size(), (pos - curr_pos) / 2));

      curr_pos = pos;
      init_support_child(pos, each_child_size[2], c_leaves[2], threshold);
      c_sizes[2] = (pos - curr_pos) / 2;

      //c.push_back(ENCODE(each_child_size[2].size(), (pos - curr_pos) / 2));

      curr_pos = pos;
      init_support_child(pos, each_child_size[3], c_leaves[3], threshold);
      c_sizes[3] = (pos - curr_pos) / 2;

      pos++;
      if(c_sizes[0] >= threshold) {
        c.push_back(child_info(each_child_size[0].size(), c_sizes[0], c_leaves[0]));
      } else {
        c.push_back(child_info(0, c_sizes[0], c_leaves[0]));
      }
      if(c_sizes[1] >= threshold) {
        c.push_back(child_info(each_child_size[1].size(), c_sizes[1], c_leaves[1]));
      } else {
        c.push_back(child_info(0, c_sizes[1], c_leaves[1]));
      }
      if(c_sizes[2] >= threshold) {
        c.push_back(child_info(each_child_size[2].size(), c_sizes[2], c_leaves[2]));
      } else {
        c.push_back(child_info(0, c_sizes[2], c_leaves[2]));
      }

      leaves = c_leaves[0] + c_leaves[1] + c_leaves[2] + c_leaves[3];

      if(c_sizes[0] >= threshold) {
        c.insert(c.end(), each_child_size[0].begin(), each_child_size[0].end());
      }
      if(c_sizes[1] >= threshold) {
        c.insert(c.end(), each_child_size[1].begin(), each_child_size[1].end());
      }
      if(c_sizes[2] >= threshold) {
        c.insert(c.end(), each_child_size[2].begin(), each_child_size[2].end());
      }
      if(c_sizes[3] >= threshold) {
        c.insert(c.end(), each_child_size[3].begin(), each_child_size[3].end());
      }
    }

    void check_size(uint64_t &pos, uint64_t &leaves, uint64_t curr_child, uint64_t curr_size, uint64_t curr_n_leaves) {
      if(pos + 3 < tree.size() && tree[pos] && tree[pos + 1] && !tree[pos + 2] && !tree[pos + 3]) {
        pos += 4;
        leaves++;
        return;
      }

      if(pos + 1 < tree.size() && tree[pos] && !tree[pos + 1]) {
        pos += 2;
        return;
      }

      if(curr_size < threshold) {
        uint64_t curr_pos = pos;
        pos++;
        check_size(pos, leaves, curr_child, curr_size, curr_n_leaves);
        check_size(pos, leaves, curr_child, curr_size, curr_n_leaves);
        check_size(pos, leaves, curr_child, curr_size, curr_n_leaves);
        check_size(pos, leaves, curr_child, curr_size, curr_n_leaves);
        pos++;

        //debug(tree_support.find_close(curr_pos) - curr_pos + 1, pos - curr_pos);
        return;
      }

      uint64_t c_size[4] = {GET_NODES(child_support[curr_child].size_tree),
                            GET_NODES(child_support[curr_child + 1].size_tree),
                            GET_NODES(child_support[curr_child + 2].size_tree),
                            0};
      uint64_t c_skip[4] = {GET_SKIPS(child_support[curr_child].size_tree),
                            GET_SKIPS(child_support[curr_child + 1].size_tree),
                            GET_SKIPS(child_support[curr_child + 2].size_tree),
                            0};
      uint64_t c_leaves[4] = {child_support[curr_child].n_leaves,
                              child_support[curr_child + 1].n_leaves,
                              child_support[curr_child + 2].n_leaves,
                              0};
      c_size[3] = curr_size - (c_size[0] + c_size[1] + c_size[2] + 1); // 3 subtrees + root
      c_leaves[3] = curr_n_leaves - (c_leaves[0] + c_leaves[1] + c_leaves[2]);

      pos++;

      uint64_t curr_pos = pos;
      uint64_t curr_leaves = leaves;
      check_size(pos, leaves, curr_child + 3, c_size[0], c_leaves[0]);
      assert(c_size[0] * 2 == pos - curr_pos);
      assert(c_leaves[0] == leaves - curr_leaves);

      curr_pos = pos;
      curr_leaves = leaves;
      check_size(pos, leaves, curr_child + 3 + c_skip[0], c_size[1], c_leaves[1]);
      assert(c_size[1] * 2 == pos - curr_pos);
      assert(c_leaves[1] == leaves - curr_leaves);

      curr_pos = pos;
      curr_leaves = leaves;
      check_size(pos, leaves, curr_child + 3 + c_skip[0] + c_skip[1], c_size[2], c_leaves[2]);
      assert(c_size[2] * 2 == pos - curr_pos);
      assert(c_leaves[2] == leaves - curr_leaves);

      curr_pos = pos;
      curr_leaves = leaves;
      check_size(pos, leaves, curr_child + 3 + c_skip[0] + c_skip[1] + c_skip[2], c_size[3], c_leaves[3]);
      assert(c_size[3] * 2 == pos - curr_pos);
      assert(c_leaves[3] == leaves - curr_leaves);

      pos++;
    }

    void traverse(uint64_t m_size, uint64_t &pos, uint64_t &size_tree, uint64_t &n_leaves) {
      assert(tree[pos] == 1);
      uint64_t depth = tree[pos]; // this should be always (

      uint64_t curr_pos = pos;
      pos++;
      while(depth != 0) {
        if(tree[pos]) {
          m_size /= 2;
          if(m_size == 1) {
            n_leaves++;
          }
        } else {
          m_size *= 2;
        }
        depth += (tree[pos] ? 1 : -1);
        pos++;
      } // always end at position +1 of last ) of this tree
      size_tree = (pos - curr_pos) / 2; 
    }

  public:
    uint64_t size() {
      uint64_t m = 0;
      sdsl::rank_support_rrr<1, 127> rank(&l);
      return rank(l.size());
    }

    uint64_t size_matrix() { return rmsize; }
    uint64_t nodes() { return tree.size() / 2; }

    k2_bp() {}
    
    k2_bp(plain_tree &pd) {
      tree = bit_vector(pd.tree.size(), 0);
      for(uint64_t i = 0; i < pd.tree.size(); i++) tree[i] = pd.tree[i];

      last_bit_t = tree.size();
      tree_support = bp_support_sada<>(&tree);

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

      uint64_t pos = 0;
      uint64_t n_leaves = 0;
      threshold = std::sqrt(tree.size() / 2);
      init_support_child(pos, child_support, n_leaves, threshold);
    }


    k2_bp(vector< pair< uint64_t, uint64_t > > &ones, uint64_t n = -1) { 
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
      //cout << "Real Size Matrix: " << rmsize << "x" << rmsize << endl;
      //cout << "Size Matrix: " << msize << "x" << msize << endl;
      //cout << "Height Tree: " << height_tree << endl;
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
      //cout << "Initialize recursion..." << endl;
#endif // DEBUG
      stack< tuple< uint64_t, uint64_t, uint64_t, uint64_t, vector< uint64_t >::iterator, uint64_t, bool, bool >,
             vector< tuple< uint64_t, uint64_t, uint64_t, uint64_t, vector< uint64_t >::iterator, uint64_t, bool, bool > > > recursion;
      recursion.push(make_tuple(msize, 0, 0, 0, ia_ones.begin(), ia_ones.size(), true, false));

      vector< uint64_t > bv_tree;
      vector< uint64_t > bv_l;

      uint64_t pos_to_add = 0;
      uint64_t pos_to_add_l = 0;

      while(!recursion.empty()) {
        auto [subm_size, init_x, init_y, smin, ia, n_ia, one_one, flag] = recursion.top();
#ifdef DEBUG
        //cout << "Recursion call..." << endl;
        //cout << "Current Sub Matrix Size: " << subm_size << " init x: " << init_x << " init y: " << init_y << " visited: " << flag << endl;
        //cout << " n ia: " << n_ia << endl;
#endif // DEBUG
        recursion.pop();
        
        if(flag) {
#ifdef DEBUG
          balance_string += ")";
          //cout << "Adding )..." << endl;
#endif // DEBUG
          add_zero(bv_tree, pos_to_add);
          continue;
        }

        recursion.push(make_tuple(subm_size, init_x, init_y, smin, ia, n_ia, one_one, true));
#ifdef DEBUG
        balance_string += "(";
        //cout << "Adding (..." << endl;
#endif // DEBUG
        add_one(bv_tree, pos_to_add);

        if(!one_one) {
          continue;
        }

        
        if(subm_size == k) {
          vector< int64_t > t(4, 0);

          for(size_t i = 0; i < n_ia; i++) {
#ifdef DEBUG
            //cout << "IA[" << i << "] = " << ia[i] << endl;
            //cout << "smin = " << smin << endl;
#endif
            int64_t pos = (int64_t) (ia[i] - smin);
#ifdef DEBUG
            //cout << "pos = " << pos << endl;
#endif
            t[pos] = 1;
          }

          for(const auto& bit : t) {
            if(bit) add_one(bv_l, pos_to_add_l);
            else add_zero(bv_l, pos_to_add_l);
          }

#ifdef DEBUG
          balance_string += "(";
          //cout << "Adding (..." << endl;
#endif // DEBUG
          add_one(bv_tree, pos_to_add);
#ifdef DEBUG
          balance_string += ")";
          //cout << "Adding )..." << endl;
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
      //cout << "Result: " << balance_string << "..." << endl;
      //cout << "Init tree " << pos_to_add << "..." << endl;
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
      //cout << "Init L " << pos_to_add_l << "..." << endl;
#endif // DEBUG
      auto aux_l = bit_vector(pos_to_add_l, 0);
      for(const auto& bit : bv_l) aux_l[bit] = 1;
      l = bv_leaves(aux_l);

#ifdef DEBUG
      //cout << "Init Tree support..." << endl;
#endif // DEBUG
      tree_support = bp_support_sada<>(&tree);

#ifdef DEBUG
      pos = 0;
      n_leaves = 0;
      check_size(pos, n_leaves, 0, tree.size() / 2, l.size() / 4);
      //cout << "End k2tree building..." << endl;
#endif // DEBUG
    }

    void add_child_info(uint64_t threshold_ = 0) {
      if(threshold_ == 0)
        threshold = std::sqrt(tree.size() / 2);
      else
        threshold = threshold_;

      uint64_t pos = 0;
      uint64_t n_leaves = 0;
      init_support_child(pos, child_support, n_leaves, threshold);
    }

    void rec_get_pos_ones(uint64_t m_size, uint64_t x, uint64_t y, uint64_t &pos, uint64_t &n_l, std::vector< pair< uint64_t, uint64_t > > &res) {
      if(pos + 1 < tree.size() && tree[pos] && !tree[pos + 1]) {
        pos += 2;
        return;
      }

      if(m_size == k) {
        pos += 4;
        if(l[n_l++])
          res.push_back({x, y});
        if(l[n_l++])
          res.push_back({x, y + 1});
        if(l[n_l++])
          res.push_back({x + 1, y});
        if(l[n_l++])
          res.push_back({x + 1, y + 1});
        return;
      }

      pos++;
      rec_get_pos_ones(m_size / 2, x, y, pos, n_l, res);
      rec_get_pos_ones(m_size / 2, x, y + m_size / 2, pos, n_l, res);
      rec_get_pos_ones(m_size / 2, x + m_size / 2, y, pos, n_l, res);
      rec_get_pos_ones(m_size / 2, x + m_size / 2, y + m_size / 2, pos, n_l, res);
      pos++;
    }

    void get_pos_ones(std::vector< pair< uint64_t, uint64_t > > &res) {
      uint64_t pos = 0;
      uint64_t n_l = 0;
      rec_get_pos_ones(msize, 0, 0, pos, n_l, res);
    }

//    vector< pair< uint64_t, uint64_t > > get_pos_ones() {
//      stack< tuple< uint8_t, uint64_t, uint64_t >, std::vector< tuple< uint8_t, uint64_t, uint64_t > > > child_visit;
//
//      uint64_t r, c;
//      r = c = 0;
//      child_visit.push({0, r, c});
//
//      uint64_t to_read_l = 0;
//      vector< pair< uint64_t, uint64_t > > ret;
//
//      for(uint64_t i = 1; i < tree.size(); i++) {
//#ifdef DEBUG
//        //cout << "Total bits " << tree.size() << endl;
//        //cout << "Reading " << i << " bit" << endl;
//#endif // DEBUG
//        if(tree[i]) {
//#ifdef DEBUG
//          //cout << "Start of subtree" << endl;
//#endif // DEBUG
//          auto [vis, r_, c_] = child_visit.top();
//          r = r_ + (vis / k) * (1 << (height_tree - child_visit.size()));
//          c = c_ + (vis % k) * (1 << (height_tree - child_visit.size()));
//          child_visit.push({0, r, c});
//#ifdef DEBUG
//          //cout << "Level " << child_visit.size() << endl;
//          //cout << "Current row: " << r << " col: " << c << endl;
//#endif // DEBUG
//        } else {
//#ifdef DEBUG
//          //cout << "End of subtree" << endl;
//#endif // DEBUG
//          if(child_visit.size() == height_tree + 1) {
//#ifdef DEBUG
//            //cout << "Last level, read real values" << endl;
//            //cout << "L size " << l.size() << " ";
//            //cout << "Current bit " << to_read_l << endl;
//#endif // DEBUG
//            for(uint64_t j = 0; j < k * k; j++) {
//              if(l[to_read_l]) {
//                auto [vis, r_, c_] = child_visit.top();
//                ret.push_back({r_ + j / k, c_ + j % k});
//              }
//              to_read_l++;
//            }
//#ifdef DEBUG
//            //cout << "Finishing reading real values" << endl;
//#endif // DEBUG
//          }
//          child_visit.pop();
//          // means we finish to read the complete tree
//          if(child_visit.size() != 0) {
//            auto [vis, r_, c_] = child_visit.top();
//            child_visit.pop();
//            child_visit.push({vis + 1, r_, c_});
//          }
//        }
//      }
//
//      return ret;
//    }

    void new_mul(k2_bp<k, bv_leaves> &b, plain_tree &c) {
      traverse_info info_a(0, 0, 0, tree.size() / 2, l.size() / 4, 1);
      traverse_info info_b(0, 0, 0, b.tree.size() / 2, b.l.size() / 4, 1);
      new_mul(msize, info_a, b, info_b, c, height_tree);
    }

    void new_mul(uint64_t m_size, traverse_info &info_a,
                 k2_bp<k, bv_leaves> &b, traverse_info &info_b,
                 plain_tree &c, uint64_t curr_h) {
      assert(tree[info_a.pos]);
      assert(b.tree[info_b.pos]);

      if(tree[info_a.pos] && !tree[info_a.pos + 1]) { // result is 0
        c.reserve(2, 0);
        c.tree.push_back(1);
        c.tree.push_back(0);
        c.height_tree = curr_h;
        c.m = m;
        c.msize = msize;
        c.rmsize = rmsize;
        info_a.pos += 2;
        info_b.pos += info_b.size << 1;
        return;
      }

      if(b.tree[info_b.pos] && !b.tree[info_b.pos + 1]) { // result is 0
        c.reserve(2, 0);
        c.tree.push_back(1);
        c.tree.push_back(0);
        c.height_tree = curr_h;
        c.m = m;
        c.msize = msize;
        c.rmsize = rmsize;
        info_b.pos += 2;
        info_a.pos += info_a.size << 1;
        return;
      }

      if(m_size == k) {
        uint8_t aux_l = minimat_mul(l.get_int(info_a.l << 2, 4), b.l.get_int(info_b.l << 2, 4));
        if(aux_l > 0) {
          c.reserve(4, 4);
          c.tree.push_back(1);
          c.tree.push_back(1);
          c.tree.push_back(0);
          c.tree.push_back(0);
          c.l.push_back(aux_l);
        } else {
          c.reserve(2, 0);
          c.tree.push_back(1);
          c.tree.push_back(0);
        }
        c.height_tree = curr_h;
        c.m = m;
        c.msize = msize;
        c.rmsize = rmsize;

        info_a.pos += 4;
        info_b.pos += 4;
        info_a.l++;
        info_b.l++;
        return;
      }

//      //  A_0 | A_1
//      //  ---------
//      //  A_2 | A_3
//      uint64_t A_0, A_1, A_2, A_3;
//      uint64_t A_0_L, A_1_L, A_2_L, A_3_L;
//
//      //  B_0 | B_1
//      //  ---------
//      //  B_2 | B_3
//      uint64_t B_0, B_1, B_2, B_3;
//      uint64_t B_0_L, B_1_L, B_2_L, B_3_L;
// 
//      //  C_0 | C_1
//      //  ---------
//      //  C_2 | C_3
//      plain_tree C_0, C_1, C_2, C_3;
//      plain_tree C_0_0, C_1_2, C_0_1, C_1_3, C_2_0, C_3_2, C_2_1, C_3_3;
      traverse_info as[4], bs[4];
      bool da, db;
      da = db = false;
      if(child_support.size() > 0 && info_a.size >= threshold) {
        as[0] = traverse_info(info_a.pos + 1, info_a.l,
                              info_a.node + 3,
                              GET_NODES(child_support[info_a.node].size_tree),
                              child_support[info_a.node].n_leaves, 1);

        as[1] = traverse_info(info_a.pos + 1 + as[0].size * 2,
                              info_a.l + as[0].n_l,
                              as[0].node + GET_SKIPS(child_support[info_a.node].size_tree),
                              GET_NODES(child_support[info_a.node + 1].size_tree),
                              child_support[info_a.node + 1].n_leaves, 1);

        as[2] = traverse_info(info_a.pos + 1 + as[0].size * 2 + as[1].size * 2,
                              info_a.l + as[0].n_l + as[1].n_l,
                              as[1].node + GET_SKIPS(child_support[info_a.node + 1].size_tree),
                              GET_NODES(child_support[info_a.node + 2].size_tree),
                              child_support[info_a.node + 2].n_leaves, 1);

        as[3] = traverse_info(info_a.pos + 1 + as[0].size * 2 + as[1].size * 2 + as[2].size * 2,
                              info_a.l + as[0].n_l + as[1].n_l + as[2].n_l,
                              as[2].node + GET_SKIPS(child_support[info_a.node + 2].size_tree),
                              info_a.size - (as[0].size + as[1].size + as[2].size + 1),
                              info_a.n_l - (as[0].n_l + as[1].n_l + as[2].n_l), 1);
      } else if(dynamic_support.size() > 0 && info_a.size >= 3) {
        as[0] = traverse_info(info_a.pos + 1, info_a.l,
                              info_a.node + 3,
                              GET_NODES(dynamic_support[info_a.node].size_tree),
                              dynamic_support[info_a.node].n_leaves, 1);

        as[1] = traverse_info(info_a.pos + 1 + as[0].size * 2,
                              info_a.l + as[0].n_l,
                              as[0].node + GET_SKIPS(dynamic_support[info_a.node].size_tree),
                              GET_NODES(dynamic_support[info_a.node + 1].size_tree),
                              dynamic_support[info_a.node + 1].n_leaves, 1);

        as[2] = traverse_info(info_a.pos + 1 + as[0].size * 2 + as[1].size * 2,
                              info_a.l + as[0].n_l + as[1].n_l,
                              as[1].node + GET_SKIPS(dynamic_support[info_a.node + 1].size_tree),
                              GET_NODES(dynamic_support[info_a.node + 2].size_tree),
                              dynamic_support[info_a.node + 2].n_leaves, 1);

        as[3] = traverse_info(info_a.pos + 1 + as[0].size * 2 + as[1].size * 2 + as[2].size * 2,
                              info_a.l + as[0].n_l + as[1].n_l + as[2].n_l,
                              as[2].node + GET_SKIPS(dynamic_support[info_a.node + 2].size_tree),
                              info_a.size - (as[0].size + as[1].size + as[2].size + 1),
                              info_a.n_l - (as[0].n_l + as[1].n_l + as[2].n_l), 1);
      } else { // traverse
        uint64_t curr_pos = info_a.pos;
        uint64_t curr_leaves = 0;
        init_support_child(curr_pos, dynamic_support, curr_leaves, 3);
        da = true;
        as[0] = traverse_info(info_a.pos + 1, info_a.l,
                              3,
                              GET_NODES(dynamic_support[0].size_tree),
                              dynamic_support[0].n_leaves, 1);

        as[1] = traverse_info(info_a.pos + 1 + as[0].size * 2,
                              info_a.l + as[0].n_l,
                              as[0].node + GET_SKIPS(dynamic_support[0].size_tree),
                              GET_NODES(dynamic_support[1].size_tree),
                              dynamic_support[1].n_leaves, 1);

        as[2] = traverse_info(info_a.pos + 1 + as[0].size * 2 + as[1].size * 2,
                              info_a.l + as[0].n_l + as[1].n_l,
                              as[1].node + GET_SKIPS(dynamic_support[1].size_tree),
                              GET_NODES(dynamic_support[2].size_tree),
                              dynamic_support[2].n_leaves, 1);

        as[3] = traverse_info(info_a.pos + 1 + as[0].size * 2 + as[1].size * 2 + as[2].size * 2,
                              info_a.l + as[0].n_l + as[1].n_l + as[2].n_l,
                              as[2].node + GET_SKIPS(dynamic_support[2].size_tree),
                              info_a.size - (as[0].size + as[1].size + as[2].size + 1),
                              info_a.n_l - (as[0].n_l + as[1].n_l + as[2].n_l), 1);
//        uint64_t curr_pos = info_a.pos + 1;
//
//        as[0].pos = info_a.pos + 1;
//        as[0].l = info_a.l;
//        traverse(m_size / 2, curr_pos, as[0].size, as[0].n_l);
//
//        as[1].pos = curr_pos;
//        as[1].l = info_a.l + as[0].n_l;
//        traverse(m_size / 2, curr_pos, as[1].size, as[1].n_l);
//
//        as[2].pos = curr_pos;
//        as[2].l = info_a.l + as[0].n_l + as[1].n_l;
//        traverse(m_size / 2, curr_pos, as[2].size, as[2].n_l);
//
//        as[3].pos = curr_pos;
//        as[3].l = info_a.l + as[0].n_l + as[1].n_l + as[2].n_l;
//        traverse(m_size / 2, curr_pos, as[3].size, as[3].n_l);
      }

      if(b.child_support.size() > 0 && info_b.size >= b.threshold) {
        bs[0] = traverse_info(info_b.pos + 1, info_b.l,
                              info_b.node + 3,
                              GET_NODES(b.child_support[info_b.node].size_tree),
                              b.child_support[info_b.node].n_leaves, 1);

        bs[1] = traverse_info(info_b.pos + 1 + bs[0].size * 2,
                              info_b.l + bs[0].n_l,
                              bs[0].node + GET_SKIPS(b.child_support[info_b.node].size_tree),
                              GET_NODES(b.child_support[info_b.node + 1].size_tree),
                              b.child_support[info_b.node + 1].n_leaves, 1);

        bs[2] = traverse_info(info_b.pos + 1 + bs[0].size * 2 + bs[1].size * 2,
                              info_b.l + bs[0].n_l + bs[1].n_l,
                              bs[1].node + GET_SKIPS(b.child_support[info_b.node + 1].size_tree),
                              GET_NODES(b.child_support[info_b.node + 2].size_tree),
                              b.child_support[info_b.node + 2].n_leaves, 1);

        bs[3] = traverse_info(info_b.pos + 1 + bs[0].size * 2 + bs[1].size * 2 + bs[2].size * 2,
                              info_b.l + bs[0].n_l + bs[1].n_l + bs[2].n_l,
                              bs[2].node + GET_SKIPS(b.child_support[info_b.node + 2].size_tree),
                              info_b.size - (bs[0].size + bs[1].size + bs[2].size + 1),
                              info_b.n_l - (bs[0].n_l + bs[1].n_l + bs[2].n_l), 1);
      } else if(b.dynamic_support.size() > 0 && info_a.size >= 3) {
        bs[0] = traverse_info(info_b.pos + 1, info_b.l,
                              info_b.node + 3,
                              GET_NODES(b.dynamic_support[info_b.node].size_tree),
                              b.dynamic_support[info_b.node].n_leaves, 1);

        bs[1] = traverse_info(info_b.pos + 1 + bs[0].size * 2,
                              info_b.l + bs[0].n_l,
                              bs[0].node + GET_SKIPS(b.dynamic_support[info_b.node].size_tree),
                              GET_NODES(b.dynamic_support[info_b.node + 1].size_tree),
                              b.dynamic_support[info_b.node + 1].n_leaves, 1);

        bs[2] = traverse_info(info_b.pos + 1 + bs[0].size * 2 + bs[1].size * 2,
                              info_b.l + bs[0].n_l + bs[1].n_l,
                              bs[1].node + GET_SKIPS(b.dynamic_support[info_b.node + 1].size_tree),
                              GET_NODES(b.dynamic_support[info_b.node + 2].size_tree),
                              b.dynamic_support[info_b.node + 2].n_leaves, 1);

        bs[3] = traverse_info(info_b.pos + 1 + bs[0].size * 2 + bs[1].size * 2 + bs[2].size * 2,
                              info_b.l + bs[0].n_l + bs[1].n_l + bs[2].n_l,
                              bs[2].node + GET_SKIPS(b.dynamic_support[info_b.node + 2].size_tree),
                              info_b.size - (bs[0].size + bs[1].size + bs[2].size + 1),
                              info_b.n_l - (bs[0].n_l + bs[1].n_l + bs[2].n_l), 1);
      } else { // traverse
        uint64_t curr_pos = info_b.pos;
        uint64_t curr_leaves = 0;
        b.init_support_child(curr_pos, b.dynamic_support, curr_leaves, 3);
        db = true;
        bs[0] = traverse_info(info_b.pos + 1, info_b.l,
                              3,
                              GET_NODES(b.dynamic_support[0].size_tree),
                              b.dynamic_support[0].n_leaves, 1);

        bs[1] = traverse_info(info_b.pos + 1 + bs[0].size * 2,
                              info_b.l + bs[0].n_l,
                              bs[0].node + GET_SKIPS(b.dynamic_support[0].size_tree),
                              GET_NODES(b.dynamic_support[1].size_tree),
                              b.dynamic_support[1].n_leaves, 1);

        bs[2] = traverse_info(info_b.pos + 1 + bs[0].size * 2 + bs[1].size * 2,
                              info_b.l + bs[0].n_l + bs[1].n_l,
                              bs[1].node + GET_SKIPS(b.dynamic_support[1].size_tree),
                              GET_NODES(b.dynamic_support[2].size_tree),
                              b.dynamic_support[2].n_leaves, 1);

        bs[3] = traverse_info(info_b.pos + 1 + bs[0].size * 2 + bs[1].size * 2 + bs[2].size * 2,
                              info_b.l + bs[0].n_l + bs[1].n_l + bs[2].n_l,
                              bs[2].node + GET_SKIPS(b.dynamic_support[2].size_tree),
                              info_b.size - (bs[0].size + bs[1].size + bs[2].size + 1),
                              info_b.n_l - (bs[0].n_l + bs[1].n_l + bs[2].n_l), 1);
//        uint64_t curr_pos = info_b.pos + 1;
//
//        bs[0].pos = info_b.pos + 1;
//        bs[0].l = info_b.l;
//        b.traverse(m_size / 2, curr_pos, bs[0].size, bs[0].n_l);
//
//        bs[1].pos = curr_pos;
//        bs[1].l = info_b.l + bs[0].n_l;
//        b.traverse(m_size / 2, curr_pos, bs[1].size, bs[1].n_l);
//
//        bs[2].pos = curr_pos;
//        bs[2].l = info_b.l + bs[0].n_l + bs[1].n_l;
//        b.traverse(m_size / 2, curr_pos, bs[2].size, bs[2].n_l);
//
//        bs[3].pos = curr_pos;
//        bs[3].l = info_b.l + bs[0].n_l + bs[1].n_l + bs[2].n_l;
//        b.traverse(m_size / 2, curr_pos, bs[3].size, bs[3].n_l);
      }

      //  C_0 | C_1
      //  ---------
      //  C_2 | C_3
      plain_tree c_[4];
      plain_tree aux_c[2];
      traverse_info save_a, save_b;
      save_a = as[0]; save_b = bs[0];
      new_mul(m_size / 2, as[0], b, bs[0], aux_c[0], curr_h - 1);
      as[0] = save_a; bs[0] = save_b;

      save_a = as[1]; save_b = bs[2];
      new_mul(m_size / 2, as[1], b, bs[2], aux_c[1], curr_h - 1);
      as[1] = save_a; bs[2] = save_b;

      aux_c[0].binsum(aux_c[1], c_[0]);
      aux_c[0].destroy();
      aux_c[1].destroy();
      aux_c[0] = aux_c[1] = plain_tree();

      save_a = as[0]; save_b = bs[1];
      new_mul(m_size / 2, as[0], b, bs[1], aux_c[0], curr_h - 1);
      as[0] = save_a; bs[1] = save_b;

      save_a = as[1]; save_b = bs[3];
      new_mul(m_size / 2, as[1], b, bs[3], aux_c[1], curr_h - 1);
      as[1] = save_a; bs[3] = save_b;

      aux_c[0].binsum(aux_c[1], c_[1]);
      aux_c[0].destroy();
      aux_c[1].destroy();
      aux_c[0] = aux_c[1] = plain_tree();

      save_a = as[2]; save_b = bs[0];
      new_mul(m_size / 2, as[2], b, bs[0], aux_c[0], curr_h - 1);
      as[2] = save_a; bs[0] = save_b;

      save_a = as[3]; save_b = bs[2];
      new_mul(m_size / 2, as[3], b, bs[2], aux_c[1], curr_h - 1);
      as[3] = save_a; bs[2] = save_b;

      aux_c[0].binsum(aux_c[1], c_[2]);
      aux_c[0].destroy();
      aux_c[1].destroy();
      aux_c[0] = aux_c[1] = plain_tree();

      save_a = as[2]; save_b = bs[1];
      new_mul(m_size / 2, as[2], b, bs[1], aux_c[0], curr_h - 1);
      as[2] = save_a; bs[1] = save_b;

      save_a = as[3]; save_b = bs[3];
      new_mul(m_size / 2, as[3], b, bs[3], aux_c[1], curr_h - 1);
      as[3] = save_a; bs[3] = save_b;

      aux_c[0].binsum(aux_c[1], c_[3]);
      aux_c[0].destroy();
      aux_c[1].destroy();
      aux_c[0] = aux_c[1] = plain_tree();


      if(da) {
        dynamic_support.clear();
      }
      if(db) {
        b.dynamic_support.clear();
      }
      if(c_[0].tree.size() == 2 &&
         c_[1].tree.size() == 2 &&
         c_[2].tree.size() == 2 &&
         c_[3].tree.size() == 2) {
        c.tree.push_back(1);
        c.tree.push_back(0);
        c.height_tree = curr_h;
        c.m = m;
        c.msize = msize;
        c.rmsize = rmsize;
        return;
      }

      c.tree.reserve(2 + c_[0].tree.size() + c_[1].tree.size() + c_[2].tree.size() + c_[3].tree.size());
      c.l.reserve(c_[0].l.size() + c_[1].l.size() + c_[2].l.size() + c_[3].l.size());
      c.tree.push_back(1);

      c.tree.concat(c_[0].tree);
      //c.l.concat(c_[0].l, 0, c_[0].l.size());
      c.l.insert(c.l.end(), c_[0].l.begin(), c_[0].l.end());
      c_[0].destroy();

      c.tree.concat(c_[1].tree);
      //c.l.concat(c_[1].l, 0, c_[1].l.size());
      c.l.insert(c.l.end(), c_[1].l.begin(), c_[1].l.end());
      c_[1].destroy();

      c.tree.concat(c_[2].tree);
      //c.l.concat(c_[2].l, 0, c_[2].l.size());
      c.l.insert(c.l.end(), c_[2].l.begin(), c_[2].l.end());
      c_[2].destroy();

      c.tree.concat(c_[3].tree);
      c.tree.push_back(0);
      //c.l.concat(c_[3].l, 0, c_[3].l.size());
      c.l.insert(c.l.end(), c_[3].l.begin(), c_[3].l.end());
      c_[3].destroy();

      c.height_tree = curr_h;
      c.m = m;
      c.msize = msize;
      c.rmsize = rmsize;


      return;

    }

    void mul(const k2_bp<k, bv_leaves> &B, plain_tree &C) {
      uint64_t A_tree, B_tree;
      A_tree = B_tree = 0;
      uint64_t A_L, B_L;
      A_L = B_L = 0;
      sdsl::int_vector<4> A_L_S(l.size() / 4, 0);
      sdsl::int_vector<4> B_L_S(B.l.size() / 4, 0);

      mul(A_tree, A_L, 0, A_L_S, B, B_tree, B_L, 0, B_L_S, C, height_tree);
    }

    void mul(uint64_t &A_tree, uint64_t &A_L, bool A_flag, sdsl::int_vector<4> &A_L_S,
             const k2_bp<k, bv_leaves> &B, uint64_t &B_tree, uint64_t &B_L, bool B_flag, sdsl::int_vector<4> &B_L_S,
             plain_tree &C,
             uint8_t curr_h) {
      // submatrix A or B full of 0's
      bool A_f0 = (!tree[A_tree + 1]);
      bool B_f0 = (!B.tree[B_tree + 1]);
      if(A_f0 && B_f0) { 
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
        uint8_t aux_l =
          minimat_mul((A_L_S[A_L >> 2] ? A_L_S[A_L >> 2] : A_L_S[A_L >> 2] = l.get_int(A_L, 4)),
                      (B_L_S[B_L >> 2] ? B_L_S[B_L >> 2] : B_L_S[B_L >> 2] = B.l.get_int(B_L, 4)));
        if(aux_l > 0) {
          C.reserve(4, 4);
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
      C.l.reserve(C_0.l.size() + C_1.l.size() + C_2.l.size() + C_3.l.size());
      C.tree.push_back(1);

      C.tree.concat(C_0.tree);
      //C.l.concat(C_0.l, 0, C_0.l.size());
      C.l.insert(C.l.end(), C_0.l.begin(), C_0.l.end());
      C_0.destroy();

      C.tree.concat(C_1.tree);
      //C.l.concat(C_1.l, 0, C_1.l.size());
      C.l.insert(C.l.end(), C_1.l.begin(), C_1.l.end());
      C_1.destroy();

      C.tree.concat(C_2.tree);
      //C.l.concat(C_2.l, 0, C_2.l.size());
      C.l.insert(C.l.end(), C_2.l.begin(), C_2.l.end());
      C_2.destroy();

      C.tree.concat(C_3.tree);
      C.tree.push_back(0);
      //C.l.concat(C_3.l, 0, C_3.l.size());
      C.l.insert(C.l.end(), C_3.l.begin(), C_3.l.end());
      C_3.destroy();

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
      uint64_t values = child_support.size();
      out.write((char*) &values, sizeof(uint64_t));
      out.write((char*) child_support.data(), values * sizeof(child_info));

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
      uint64_t size;
      in.read((char*) &size, sizeof(uint64_t));
      child_support.resize(size, child_info());
      in.read((char*) child_support.data(), size * sizeof(child_info));

      sdsl::load(leaves, in);
      //leaves.load(in);
      rank_leaves.load(in, &leaves);

      tree.load(in);
      tree_support.load(in, &tree);

      //l.load(in);
      sdsl::load(l, in);
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
      cout << "  Tree        : " << (size_in_bytes(tree)) * 8 << "," << (double) (size_in_bytes(tree)) * 8 / size() << "," << (double) (size_in_bytes(tree)) * 8 / total << endl;
      cout << "  Tree Support: " << (size_in_bytes(tree_support)) * 8 << "," << (double) (size_in_bytes(tree_support)) * 8 / size() << "," << (double) (size_in_bytes(tree_support)) * 8 / total << endl;
      cout << "  L           : " << (size_in_bytes(l)) * 8 << "," << (double) (size_in_bytes(l)) * 8 / size() << "," << (double) (size_in_bytes(l)) * 8 / total << endl;
      cout << "  child supp  : " << child_support.size() * sizeof(child_info) * 8 << "," << (double) (child_support.size() * sizeof(child_info) * 8) / size() << "," << (double) child_support.size() * sizeof(child_info) * 8 / total << endl;
      cout << "  leaves      : " << (size_in_bytes(leaves) + size_in_bytes(rank_leaves)) * 8 << "," << (double) (size_in_bytes(leaves) + size_in_bytes(rank_leaves)) * 8 / size() << "," << (double) (size_in_bytes(leaves) + size_in_bytes(rank_leaves)) * 8 / total << endl;
#endif
      return total;
    }

    friend ostream& operator<<(ostream& os, const k2_bp<k, bv_leaves> &k2tree) {
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
