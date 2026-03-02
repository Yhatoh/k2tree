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
template< uint64_t k = 2, class bv_leaves = bit_vector, uint64_t bexc = 256 >
class k2_bp {
  public:
    uint64_t height_tree;

    sdsl::int_vector<> exc_min_samples;
    sdsl::int_vector<> exc_samples;
    sdsl::int_vector<> leaves_samples;

    std::vector< child_info > child_support;
    std::vector< child_info > dynamic_support;
    bit_vector tree; // k2tree
    uint64_t last_bit_t; // universe

    bv_leaves l; // real values
    uint64_t last_bit_l; // universe

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
      if(pos + 3 < tree.size() && tree.get_int(pos, 4) == 3) {
        pos += 4;
        leaves++;
        return;
      }

      if(pos + 1 < tree.size() && tree.get_int(pos, 2) == 1) {
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

    inline uint64_t count(uint64_t num) {
      uint64_t x = num;
      uint64_t y = num >> 1;
      uint64_t hi = x & y;
      uint64_t lo = ~ (x | y);

      uint64_t bits = (hi & (lo >> 2)) & 2305843009213693951;
      return __builtin_popcountll(bits);
    }

    inline uint64_t rank_leave(uint64_t i, uint64_t j) {

      uint64_t ret = 0;
      uint64_t bit = i;
      for(; bit + 64 < j; bit += 64) {
        uint64_t extra = 0;
        uint64_t read = tree.get_int(bit, 64);
        if(bit + 64 < j) {
          uint64_t len = (6 > tree.size() - (bit + 61) ?
                          tree.size() - (bit + 61) :
                          6);
          extra = count(tree.get_int(bit + 61, len) | (((uint64_t) -1) << len));
        }
        ret += count(read) + extra;
      }
      if(j > bit) {
        uint64_t dist = j - bit;
        if(dist == 64) ret += count(tree.get_int(bit, 64));
        else ret += count(tree.get_int(bit, dist) | (((uint64_t) -1) << (dist)));

        uint64_t extra = 0;

        if(dist == 1) {
          uint64_t len_extra = (j - 1 + 4 <= tree.size() ? 4 : tree.size() - (j - 1));
          uint64_t read = tree.get_int(j - 1, len_extra);
          extra = 
            count(read | (((uint64_t) -1) << len_extra));
          ret += extra;
        } else if(dist == 2) {
          uint64_t len_extra = (j - 2 + 5 <= tree.size() ? 5 : tree.size() - (j - 2));
          uint64_t read = tree.get_int(j - 2, len_extra);
          extra = 
            count(read | (((uint64_t) -1) << len_extra));
          ret += extra;
        } else {
          uint64_t len_extra = (j + 3 <= tree.size() ? 6 : tree.size() - (j - 3));
          uint64_t read = tree.get_int(j - 3, len_extra);
          extra = 
            count(read | (((uint64_t) -1) << len_extra));
          ret += extra;
        }
      }
      return ret;
    }

    inline void fasttraverse(uint64_t &pos, uint64_t &size_tree, uint64_t &n_leaves) {
      if(tree.get_int(pos, 4) == 3) {
        pos += 4; size_tree = 2; n_leaves = 1;
        return;
      }
      if(tree.get_int(pos, 2) == 1) {
        pos += 2; size_tree = 1; n_leaves = 0;
        return;
      }

      uint64_t curr_pos = pos;
      int64_t depth = 0;

      uint64_t bits = tree.get_int(pos, 11); // 12 is the smaller size it should receive
      n_leaves += count(bits | ((uint64_t)-1 << 11));

      int64_t plus_one = __builtin_popcountll(bits & ((1ULL << 8) - 1ULL));
      depth = plus_one - (8 - plus_one);

      pos += 8;

      while(depth > 0) {
        bits = tree.get_int(pos, 11);
        uint8_t end__= end_tree[depth - 1][bits & ((1LL << 8) - 1LL)];
        if(end__== 8) {
          plus_one = __builtin_popcountll(bits & ((1LL << 8) - 1LL));
          depth += plus_one - (8 - plus_one);
          n_leaves += count(bits | ((uint64_t) -1 << 11));
          pos += 8;
        } else {
          //plus_one = __builtin_popcountll(bits & ((1LL << (end__+ 1)) - 1LL));
          depth = 0;
          n_leaves += count(bits | ((uint64_t) -1 << (end__+ 1)));
          pos += end__+ 1;
        }
      }
      size_tree = (pos - curr_pos) / 2;
    }

    inline void ultratraverse(uint64_t &pos, int64_t excess, uint64_t &size_tree, uint64_t &n_leaves) {
      if(tree.get_int(pos, 4) == 3) {
        pos += 4; size_tree = 2; n_leaves = 1;
        return;
      }
      if(tree.get_int(pos, 2) == 1) {
        pos += 2; size_tree = 1; n_leaves = 0;
        return;
      }

      uint64_t curr_pos = pos;
      int64_t obj_excess = excess;
      pos++;
      uint64_t obj = bexc - pos % bexc + pos;
      for(; pos + 16 < obj; pos += 16) {
        // found micro block
        uint64_t bits = tree.get_int(pos, 16);
        if(obj_excess >= excess + exc_min_micro[bits] + 1) { 
          for(uint8_t i = 0; i < 16; bits = bits >> 1, i++) { // reading bit
            pos++;
            if(bits & 1) excess++;
            else excess--;
           
            if(obj_excess == excess + 1) {
              size_tree = (pos - curr_pos) / 2;
              n_leaves = rank_leave(curr_pos, pos);
              return;
            }
          }
        }

        excess += exc_micro[bits];
      }

      // read last part;
      {
        uint64_t extra = bexc - pos % bexc;
        uint64_t bits = tree.get_int(pos, extra);
        for(uint8_t i = 0; i < extra; bits = bits >> 1, i++) { // reading bit
          pos++;
          if(bits & 1) excess++;
          else excess--;

          if(obj_excess == excess + 1) {
            size_tree = (pos - curr_pos) / 2;
            n_leaves = rank_leave(curr_pos, pos);
            return;
          }
        }
      }

      uint64_t block = pos / bexc;
      for(;; pos += bexc) {
        // found the block
        if(excess > exc_min_samples[block]) {
          for(;; pos += 16) {
            uint64_t bits = tree.get_int(pos, 16);
            // found microblock
            if(obj_excess >= excess + exc_min_micro[bits] + 1) {
              for(uint8_t i = 0; i < 16; bits = bits >> 1, i++) {
                pos++;
                if(bits & 1) excess++;
                else excess--;

                if(obj_excess == excess + 1) {
                  size_tree = (pos - curr_pos) / 2;
                  n_leaves = rank_leave(curr_pos, pos);
                  return;
                }
              }
            }
            excess += exc_micro[bits];
          }
        }
        excess = exc_samples[block++];
      }
    }

    void build_exc_sample() {

      int64_t excess = 1;
      int64_t min_excess = 1;
      exc_min_samples.resize((tree.size() + bexc - 1) / bexc + 1);
      exc_samples.resize((tree.size() + bexc - 1) / bexc + 1);
      uint64_t block = 0;
      for(size_t i = 1; i < tree.size(); i++) {
        if(i % bexc == 0) {
          exc_min_samples[block] = min_excess;
          exc_samples[block] = excess;
          block++;
          min_excess = LLONG_MAX;
        }
        excess += (tree[i] ? 1 : -1);
        if(excess < min_excess) min_excess = excess;
      }
      exc_min_samples[block] = min_excess;
      exc_samples[block] = excess;
      sdsl::util::bit_compress(exc_min_samples);
      sdsl::util::bit_compress(exc_samples);
    }

  public:
    uint64_t size() { return m; }

    uint64_t size_matrix() { return rmsize; }
    uint64_t nodes() { return tree.size() / 2; }

    k2_bp() {threshold = 0;}
    
    k2_bp(plain_tree &pd) {
      tree = bit_vector(pd.tree.size(), 0);
      for(uint64_t i = 0; i < pd.tree.size(); i++) tree[i] = pd.tree[i];

      m = 0;
      bit_vector aux_l = bit_vector(pd.l.size() * 4, 0);
      for(uint64_t i = 0; i < pd.l.size(); i++) {
        for(uint64_t j = 0; j < 4; j++) {
          if(pd.l[i] & (1 << j)) {
            aux_l[i * 4 + j] = 1;
            m++;
          }
        }
      }

      l = bv_leaves(aux_l);
      last_bit_l = l.size();

      height_tree = pd.height_tree;
      msize = pd.msize;
      rmsize = pd.rmsize;
      build_exc_sample();
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
      
      vector< uint64_t > ia_ones;
      for(const auto& one : ones) {
        ia_ones.push_back(bits_interleave(one.first, one.second));
      }

      sort(ia_ones.begin(), ia_ones.end());
      stack< tuple< uint64_t, uint64_t, uint64_t, uint64_t, vector< uint64_t >::iterator, uint64_t, bool, bool >,
             vector< tuple< uint64_t, uint64_t, uint64_t, uint64_t, vector< uint64_t >::iterator, uint64_t, bool, bool > > > recursion;
      recursion.push(make_tuple(msize, 0, 0, 0, ia_ones.begin(), ia_ones.size(), true, false));

      vector< uint64_t > bv_tree;
      vector< uint64_t > bv_l;

      uint64_t pos_to_add = 0;
      uint64_t pos_to_add_l = 0;

      while(!recursion.empty()) {
        auto [subm_size, init_x, init_y, smin, ia, n_ia, one_one, flag] = recursion.top();
        recursion.pop();
        
        if(flag) {
          add_zero(bv_tree, pos_to_add);
          continue;
        }

        recursion.push(make_tuple(subm_size, init_x, init_y, smin, ia, n_ia, one_one, true));
        add_one(bv_tree, pos_to_add);

        if(!one_one) {
          continue;
        }

        
        if(subm_size == k) {
          vector< int64_t > t(4, 0);

          for(size_t i = 0; i < n_ia; i++) {
            int64_t pos = (int64_t) (ia[i] - smin);
            t[pos] = 1;
          }

          for(const auto& bit : t) {
            if(bit) add_one(bv_l, pos_to_add_l);
            else add_zero(bv_l, pos_to_add_l);
          }

          add_one(bv_tree, pos_to_add);
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

      tree = bit_vector(pos_to_add, 0);
      for(const auto& bit : bv_tree) tree[bit] = 1;

      auto aux_l = bit_vector(pos_to_add_l, 0);
      for(const auto& bit : bv_l) aux_l[bit] = 1;
      l = bv_leaves(aux_l);
      
      build_exc_sample();
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

    void sum(k2_bp<k, bv_leaves> &b, k2_bp<k, bv_leaves> &c) {
      traverse_info info_a(0, 0, 0, tree.size() / 2, l.size() / 4, 1);
      traverse_info info_b(0, 0, 0, b.tree.size() / 2, b.l.size() / 4, 1);
      c.msize = msize;
      c.height_tree = height_tree;
      c.rmsize = rmsize;
      sum(msize, info_a, b, info_b, c, height_tree, 1);
    }

    void sum(uint64_t m_size, traverse_info &info_a,
             k2_bp<k, bv_leaves> &b, traverse_info &info_b,
             k2_bp<k, bv_leaves> &c, uint64_t curr_h, int64_t excess) {
      assert(tree[info_a.pos]);
      assert(b.tree[info_b.pos]);
      if(tree.get_int(info_a.pos, 2) == 1) {// copy subtree of b
        if(info_b.size == 0) {
          uint64_t curr_pos = info_b.pos;
          b.ultratraverse(curr_pos, excess, info_b.size, info_b.n_l);
        }
        c.tree.insert(c.tree.end(), b.tree.begin() + info_b.pos, b.tree.begin() + info_b.pos + (info_b.size << 1));
        c.l.insert(c.l.end(), b.l.begin() + (info_b.l << 2), b.l.begin() + ((info_b.l + info_b.n_l) << 2));
        info_a.pos += 2;
        info_b.pos += info_b.size << 1;
        info_b.l += info_b.n_l;
        return;
      }

      if(b.tree.get_int(info_b.pos, 2) == 1) {// copy subtree of a
        if(info_a.size == 0) {
          uint64_t curr_pos = info_a.pos;
          ultratraverse(curr_pos, excess, info_a.size, info_a.n_l);
        }
        c.tree.insert(c.tree.end(), tree.begin() + info_a.pos, tree.begin() + info_a.pos + (info_a.size << 1));
        c.l.insert(c.l.end(), l.begin() + (info_a.l << 2), l.begin() + ((info_a.l + info_a.n_l) << 2));
        info_b.pos += 2;
        info_a.pos += info_a.size << 1;
        info_a.l += info_a.n_l;
        return;
      }

      if(tree.get_int(info_a.pos, 4) == 3 && b.tree.get_int(info_b.pos, 4) == 3) {
        c.tree.insert(c.tree.end(), tree.begin() + info_a.pos, tree.begin() + info_a.pos + 4);
        c.l.push_back(l[(info_a.l << 2)] | b.l[(info_b.l << 2)]);
        c.l.push_back(l[(info_a.l << 2) + 1] | b.l[(info_b.l << 2) + 1]);
        c.l.push_back(l[(info_a.l << 2) + 2] | b.l[(info_b.l << 2) + 2]);
        c.l.push_back(l[(info_a.l << 2) + 3] | b.l[(info_b.l << 2) + 3]);
        info_b.pos += 4;
        info_a.pos += 4;
        info_b.l += 1;
        info_a.l += 1;
        return;
      }

      c.tree.push_back(1);
      // first submatrix
      uint64_t accum_size_a = 0;
      uint64_t accum_size_b = 0;
      uint64_t accum_l_a = 0;
      uint64_t accum_l_b = 0;

      traverse_info aux_a(info_a.pos + 1, info_a.l, 0, 0, 0, 0);
      if(info_a.size >= threshold) {
        aux_a.size = GET_NODES(child_support[info_a.node].size_tree);
        aux_a.node = info_a.node + 3;
        aux_a.n_l = child_support[info_a.node].n_leaves;
      }
      traverse_info aux_b(info_b.pos + 1, info_b.l, 0, 0, 0, 0);
      if(info_b.size >= b.threshold) {
        aux_b.size = GET_NODES(b.child_support[info_b.node].size_tree);
        aux_b.node = info_b.node + 3;
        aux_b.n_l = b.child_support[info_b.node].n_leaves;
      }
      sum(m_size / 2, aux_a, b, aux_b, c, curr_h - 1, excess + 1);
      accum_size_a += aux_a.pos - (info_a.pos + 1);
      accum_size_b += aux_b.pos - (info_b.pos + 1);
      accum_l_a += aux_a.l - info_a.l;
      accum_l_b += aux_b.l - info_b.l;
      info_a.pos = aux_a.pos;
      info_b.pos = aux_b.pos;
      info_a.l = aux_a.l;
      info_b.l = aux_b.l;

      // second submatrix
      aux_a = traverse_info(aux_a.pos, aux_a.l, 0, 0, 0, 0);
      if(info_a.size >= threshold) {
        aux_a.node = info_a.node + 3 + GET_SKIPS(child_support[info_a.node].size_tree);
        aux_a.size = GET_NODES(child_support[info_a.node + 1].size_tree);
        aux_a.n_l = child_support[info_a.node + 1].n_leaves;
      }
      aux_b = traverse_info(aux_b.pos, aux_b.l, 0, 0, 0, 0);
      if(info_b.size >= b.threshold) {
        aux_b.size = GET_NODES(b.child_support[info_b.node + 1].size_tree);
        aux_b.node = info_b.node + 3 + GET_SKIPS(b.child_support[info_b.node].size_tree);
        aux_b.n_l = b.child_support[info_b.node + 1].n_leaves;
      }
      sum(m_size / 2, aux_a, b, aux_b, c, curr_h - 1, excess + 1);
      accum_size_a += aux_a.pos - info_a.pos;
      accum_size_b += aux_b.pos - info_b.pos;
      accum_l_a += aux_a.l - info_a.l;
      accum_l_b += aux_b.l - info_b.l;
      info_a.pos = aux_a.pos;
      info_b.pos = aux_b.pos;
      info_a.l = aux_a.l;
      info_b.l = aux_b.l;

      // third submatrix
      aux_a = traverse_info(aux_a.pos, aux_a.l, 0, 0, 0, 0);
      if(info_a.size >= threshold) {
        aux_a.node = info_a.node + 3 + GET_SKIPS(child_support[info_a.node].size_tree)
                                     + GET_SKIPS(child_support[info_a.node + 1].size_tree);
        aux_a.size = GET_NODES(child_support[info_a.node + 2].size_tree);
        aux_a.n_l = child_support[info_a.node + 2].n_leaves;
      }
      aux_b = traverse_info(aux_b.pos, aux_b.l, 0, 0, 0, 0);
      if(info_b.size >= b.threshold) {
        aux_b.size = GET_NODES(b.child_support[info_b.node + 2].size_tree);
        aux_b.node = info_b.node + 3 + GET_SKIPS(b.child_support[info_b.node].size_tree)
                                     + GET_SKIPS(b.child_support[info_b.node + 1].size_tree);
        aux_b.n_l = b.child_support[info_b.node + 2].n_leaves;
      }
      sum(m_size / 2, aux_a, b, aux_b, c, curr_h - 1, excess + 1);
      accum_size_a += aux_a.pos - info_a.pos;
      accum_size_b += aux_b.pos - info_b.pos;
      accum_l_a += aux_a.l - info_a.l;
      accum_l_b += aux_b.l - info_b.l;
      info_a.pos = aux_a.pos;
      info_b.pos = aux_b.pos;
      info_a.l = aux_a.l;
      info_b.l = aux_b.l;

      // fourth submatrix
      aux_a = traverse_info(aux_a.pos, aux_a.l, 0, 0, 0, 0);
      if(info_a.size >= threshold) {
        aux_a.node = info_a.node + 3 + GET_SKIPS(child_support[info_a.node].size_tree)
                                     + GET_SKIPS(child_support[info_a.node + 1].size_tree)
                                     + GET_SKIPS(child_support[info_a.node + 2].size_tree);
        aux_a.size = info_a.size - accum_size_a / 2 - 1;
        aux_a.n_l = info_a.n_l - accum_l_a;
      }
      aux_b = traverse_info(aux_b.pos, aux_b.l, 0, 0, 0, 0);
      if(info_b.size >= b.threshold) {
        aux_b.size = info_b.size - accum_size_b / 2 - 1;
        aux_b.node = info_b.node + 3 + GET_SKIPS(b.child_support[info_b.node].size_tree)
                                     + GET_SKIPS(b.child_support[info_b.node + 1].size_tree)
                                     + GET_SKIPS(b.child_support[info_b.node + 2].size_tree);
        aux_b.n_l = info_b.n_l - accum_l_b;
      }
      sum(m_size / 2, aux_a, b, aux_b, c, curr_h - 1, excess + 1);
      c.tree.push_back(0);
      accum_size_a += aux_a.pos - info_a.pos;
      accum_size_b += aux_b.pos - info_b.pos;
      accum_l_a += aux_a.l - info_a.l;
      accum_l_b += aux_b.l - info_b.l;
      info_a.pos = aux_a.pos;
      info_b.pos = aux_b.pos;
      info_a.l = aux_a.l;
      info_b.l = aux_b.l;

      info_b.pos = aux_b.pos + 1;
      info_a.pos = aux_a.pos + 1;
    }

    void mul(k2_bp<k, bv_leaves> &b, k2_bp<k, bv_leaves> &c) {
      traverse_info info_a(0, 0, 0, tree.size() / 2, l.size() / 4, 1);
      traverse_info info_b(0, 0, 0, b.tree.size() / 2, b.l.size() / 4, 1);
      //dynamic_support.reserve(threshold * 3);
      //b.dynamic_support.reserve(b.threshold * 3);
      mul(msize, info_a, b, info_b, c, height_tree, 1);
    }

    void mul(uint64_t m_size, traverse_info &info_a,
                 k2_bp<k, bv_leaves> &b, traverse_info &info_b,
                 k2_bp<k, bv_leaves> &c, uint64_t curr_h, int64_t excess) {
      assert(tree[info_a.pos]);
      assert(b.tree[info_b.pos]);

      if(tree.get_int(info_a.pos, 2) == 1) { // result is 0
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

      if(b.tree.get_int(info_b.pos, 2) == 1) { // result is 0
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
        uint8_t aux_l = table_mul[l.get_int(info_a.l << 2, 4)][b.l.get_int(info_b.l << 2, 4)];
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
      } else {
        uint64_t curr_pos = info_a.pos + 1;

        as[0].pos = info_a.pos + 1;
        as[0].l = info_a.l;
        ultratraverse(curr_pos, excess + 1, as[0].size, as[0].n_l);

        as[1].pos = curr_pos;
        as[1].l = info_a.l + as[0].n_l;

        as[2].pos = curr_pos;
        as[2].l = info_a.l + as[0].n_l + as[1].n_l;
        ultratraverse(curr_pos, excess + 1, as[2].size, as[2].n_l);
        
        as[3].pos = curr_pos;
        as[3].l = info_a.l + as[0].n_l + as[1].n_l + as[2].n_l;
        ultratraverse(curr_pos, excess + 1, as[3].size, as[3].n_l);
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
      } else { // traverse
        uint64_t curr_pos = info_b.pos + 1;

        bs[0].pos = info_b.pos + 1;
        bs[0].l = info_b.l;
        b.ultratraverse(curr_pos, excess + 1, bs[0].size, bs[0].n_l);

        bs[1].pos = curr_pos;
        bs[1].l = info_b.l + bs[0].n_l;
        b.ultratraverse(curr_pos, excess + 1, bs[1].size, bs[1].n_l);

        bs[2].pos = curr_pos;
        bs[2].l = info_b.l + bs[0].n_l + bs[1].n_l;
        b.ultratraverse(curr_pos, excess + 1, bs[2].size, bs[2].n_l);

        bs[3].pos = curr_pos;
        bs[3].l = info_b.l + bs[0].n_l + bs[1].n_l + bs[2].n_l;
        b.ultratraverse(curr_pos, excess + 1, bs[3].size, bs[3].n_l);
      }

      //  C_0 | C_1
      //  ---------
      //  C_2 | C_3
      k2_bp<k, bv_leaves> c_[4];
      k2_bp<k, bv_leaves> aux_c[2];

      traverse_info save_a, save_b;
      save_a = as[0]; save_b = bs[0];
      mul(m_size / 2, as[0], b, bs[0], aux_c[0], curr_h - 1, excess + 1);
      as[0] = save_a; bs[0] = save_b;

      save_a = as[1]; save_b = bs[2];
      mul(m_size / 2, as[1], b, bs[2], aux_c[1], curr_h - 1, excess + 1);
      as[1] = save_a; bs[2] = save_b;

      c_[0].tree.reserve(aux_c[0].tree.size() + aux_c[1].tree.size());
      c_[0].l.reserve(aux_c[0].l.size() + aux_c[1].l.size());
      aux_c[0].sum(aux_c[1], c_[0]);
      //aux_c[0].destroy();
      //aux_c[1].destroy();
      aux_c[0] = aux_c[1] = k2_bp<k, bv_leaves>();

      save_a = as[0]; save_b = bs[1];
      mul(m_size / 2, as[0], b, bs[1], aux_c[0], curr_h - 1, excess + 1);
      as[0] = save_a; bs[1] = save_b;

      save_a = as[1]; save_b = bs[3];
      mul(m_size / 2, as[1], b, bs[3], aux_c[1], curr_h - 1, excess + 1);
      as[1] = save_a; bs[3] = save_b;

      c_[1].tree.reserve(aux_c[0].tree.size() + aux_c[1].tree.size());
      c_[1].l.reserve(aux_c[0].l.size() + aux_c[1].l.size());
      aux_c[0].sum(aux_c[1], c_[1]);
      //aux_c[0].destroy();
      //aux_c[1].destroy();
      aux_c[0] = aux_c[1] = k2_bp<k, bv_leaves>();

      save_a = as[2]; save_b = bs[0];
      mul(m_size / 2, as[2], b, bs[0], aux_c[0], curr_h - 1, excess + 1);
      as[2] = save_a; bs[0] = save_b;

      save_a = as[3]; save_b = bs[2];
      mul(m_size / 2, as[3], b, bs[2], aux_c[1], curr_h - 1, excess + 1);
      as[3] = save_a; bs[2] = save_b;

      c_[2].tree.reserve(aux_c[0].tree.size() + aux_c[1].tree.size());
      c_[2].l.reserve(aux_c[0].l.size() + aux_c[1].l.size());
      aux_c[0].sum(aux_c[1], c_[2]);
      //aux_c[0].destroy();
      //aux_c[1].destroy();
      aux_c[0] = aux_c[1] = k2_bp<k, bv_leaves>();

      save_a = as[2]; save_b = bs[1];
      mul(m_size / 2, as[2], b, bs[1], aux_c[0], curr_h - 1, excess + 1);
      as[2] = save_a; bs[1] = save_b;

      save_a = as[3]; save_b = bs[3];
      mul(m_size / 2, as[3], b, bs[3], aux_c[1], curr_h - 1, excess + 1);
      as[3] = save_a; bs[3] = save_b;


      c_[3].tree.reserve(aux_c[0].tree.size() + aux_c[1].tree.size());
      c_[3].l.reserve(aux_c[0].l.size() + aux_c[1].l.size());
      aux_c[0].sum(aux_c[1], c_[3]);
      //aux_c[0].destroy();
      //aux_c[1].destroy();
      aux_c[0] = aux_c[1] = k2_bp<k, bv_leaves>();

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

      c.tree.insert(c.tree.begin(), c_[0].tree.begin(), c_[0].tree.end());
      c.l.insert(c.l.end(), c_[0].l.begin(), c_[0].l.end());
      //c_[0].destroy();
      c_[0] = k2_bp<k, bv_leaves>();

      c.tree.insert(c.tree.begin(), c_[1].tree.begin(), c_[1].tree.end());
      c.l.insert(c.l.end(), c_[1].l.begin(), c_[1].l.end());
      //c_[1].destroy();
      c_[1] = k2_bp<k, bv_leaves>();

      c.tree.insert(c.tree.begin(), c_[2].tree.begin(), c_[2].tree.end());
      c.l.insert(c.l.end(), c_[2].l.begin(), c_[2].l.end());
      //c_[2].destroy();
      c_[2] = k2_bp<k, bv_leaves>();

      c.tree.insert(c.tree.begin(), c_[3].tree.begin(), c_[3].tree.end());
      c.tree.push_back(0);
      c.l.insert(c.l.end(), c_[3].l.begin(), c_[3].l.end());
      //c_[3].destroy();
      c_[3] = k2_bp<k, bv_leaves>();

      c.height_tree = curr_h;
      c.m = m;
      c.msize = msize;
      c.rmsize = rmsize;


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
      out.write((char*) &threshold, sizeof(uint64_t));
      uint64_t values = child_support.size();
      out.write((char*) &values, sizeof(uint64_t));
      out.write((char*) child_support.data(), values * sizeof(child_info));

      tree.serialize(out);
      l.serialize(out);
      exc_min_samples.serialize(out);
      exc_samples.serialize(out);
    }

    void load(ifstream& in) {
      // writing integers first
      in.read((char*) &msize, sizeof(uint64_t));
      in.read((char*) &rmsize, sizeof(uint64_t));
      in.read((char*) &m, sizeof(uint64_t));
      in.read((char*) &height_tree, sizeof(uint64_t));
      in.read((char*) &last_bit_t, sizeof(uint64_t));
      in.read((char*) &last_bit_l, sizeof(uint64_t));
      in.read((char*) &threshold, sizeof(uint64_t));
      uint64_t size;
      in.read((char*) &size, sizeof(uint64_t));
      child_support.resize(size, child_info());
      in.read((char*) child_support.data(), size * sizeof(child_info));

      tree.load(in);

      sdsl::load(l, in);
      sdsl::load(exc_min_samples,in);
      sdsl::load(exc_samples,in);
    }

    uint64_t size_in_bits() {
      uint64_t total = sizeof(uint64_t) * 5 +
             size_in_bytes(tree) * 8 +
             child_support.size() * sizeof(child_info) * 8 +
             size_in_bytes(l) * 8 + 
             size_in_bytes(exc_min_samples) * 8 + size_in_bytes(exc_samples) * 8 +
             65536 * 8 * 2;
#ifdef INFO_SPACE
      vector< child_info > aux_leaves; uint64_t pos = 0; uint64_t leaves = 0;
      init_support_child(pos, aux_leaves, leaves, nodes() / 2);
      cout << "Leaves:" << leaves << endl;
      cout << "BITS" << endl;
      cout << "  Tree        : " << (size_in_bytes(tree)) * 8 << "," << (double) (size_in_bytes(tree)) * 8 / size() << "," << (double) (size_in_bytes(tree)) * 8 / total << endl;
      cout << "  L           : " << (size_in_bytes(l)) * 8 << "," << (double) (size_in_bytes(l)) * 8 / size() << "," << (double) (size_in_bytes(l)) * 8 / total << endl;
      cout << "  child supp  : " << child_support.size() * sizeof(child_info) * 8 << "," << (double) (child_support.size() * sizeof(child_info) * 8) / size() << "," << (double) child_support.size() * sizeof(child_info) * 8 / total << endl;
      cout << "  microtable  : " << 65536*8*2 << "," << (double) (65536*8*2) / size() << "," << (double) 65536*8*2 / total << endl;
      cout << "  exc samples : " << size_in_bytes(exc_min_samples) * 8 + size_in_bytes(exc_samples) * 8 << "," << (double) (size_in_bytes(exc_min_samples) * 8 + size_in_bytes(exc_samples) * 8) / size() << "," << (double) (size_in_bytes(exc_min_samples) * 8 + size_in_bytes(exc_samples) * 8) / total << std::endl;

#endif
      return total;
    }
};
#endif // !K2_TREE_BP_SDSL
