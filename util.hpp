#ifndef UTIL_HPP
#define UTIL_HPP

#include <cstdint>
#include <iostream> 
#include <queue> 
#include <vector>
#include <algorithm> 
#include <utility>       
#include <map>
#include <cmath>
#include <inttypes.h>
#include <fstream>
#include <string>
#include <tuple>

#include <sdsl/vectors.hpp>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/coder_elias_delta.hpp>
#include <sdsl/suffix_arrays.hpp>

#include "block_element.hpp"

using namespace std;
// from succint repository
static const uint8_t debruijn64_mapping[64] = {
  63,  0, 58,  1, 59, 47, 53,  2,
  60, 39, 48, 27, 54, 33, 42,  3,
  61, 51, 37, 40, 49, 18, 28, 20,
  55, 30, 34, 11, 43, 14, 22,  4,
  62, 57, 46, 52, 38, 26, 32, 41,
  50, 36, 17, 19, 29, 10, 13, 21,
  56, 45, 25, 31, 35, 16,  9, 12,
  44, 24, 15,  8, 23,  7,  6,  5
};

static const uint64_t debruijn64 = 0x07EDD5E59A4E28C2ULL;

inline uint8_t bit_position(uint64_t x){
    return debruijn64_mapping[(x * debruijn64) >> 58];
}

inline uint8_t msb(uint64_t x, unsigned long& ret){
  if (!x)
    return false;

  x |= x >> 1;
  x |= x >> 2;
  x |= x >> 4;
  x |= x >> 8;
  x |= x >> 16;
  x |= x >> 32;

  x ^= x >> 1;
  ret = bit_position(x);

  return true;
}

inline uint8_t msb(uint64_t x){
  unsigned long ret = -1U;
  msb(x, ret);
  return (uint8_t)ret;
}

inline uint64_t ceil_log2(uint64_t x) {
  return (x > 1) ? msb(x - 1) + 1 : 0;
}

inline uint64_t floor_log2(uint64_t x) {
  return (x > 1) ? msb(x) : 0;
}


typedef std::pair<float, uint64_t>   heap_node; 
typedef std::tuple<uint64_t, uint64_t, std::string>   heap_node_str; 

typedef std::pair<uint64_t, int64_t> tree_node;
typedef std::tuple<uint64_t, int64_t, std::string> tree_node_str;

template< uint16_t w = 16 >
class tunstall_coder {
public:
    std::vector<std::vector<uint64_t>> D;  // Dictionary
    std::vector<uint32_t> map_table;
    //std::vector<uint16_t> compressed_seq;
    sdsl::int_vector<w> compressed_seq;
    std::vector<blockElement> block;
    //enc_vector<> block_info_prefix_sum;
    //enc_vector<> block_info_starting_position;
    uint32_t bSize;
    uint32_t sigma;
    //uint64_t D_size=65536 by default
    uint64_t D_size;

    void traverse(std::vector<tree_node>& tree, uint64_t curnode, uint64_t& curindex, uint64_t sigma,
        std::vector<uint64_t>& currcode) {
      uint64_t i, cursum;

      for (i = 0; i < sigma; ++i, ++curnode) {
        currcode.push_back(i);

        if (tree[curnode].second == -1) {
          tree[curnode].first = curindex;  // dictionary index of this code
          cursum = 0;
          for (uint64_t j = 0; j < currcode.size(); ++j) {
            cursum += map_table[currcode[j]];
            D[curindex].push_back(cursum);
          }
          curindex++;
        } else {
          tree[curnode].first = curindex;
          cursum = 0;
          for (uint64_t j = 0; j < currcode.size(); ++j) {
            cursum += map_table[currcode[j]];
            D[curindex].push_back(cursum);
          }
          curindex++;
          traverse(tree, tree[curnode].second, curindex, sigma, currcode);  // goes to children of curnode
        }
        currcode.pop_back();
      }

    }

    void traverse(std::vector<tree_node_str>& tree, uint64_t curnode, uint64_t& curindex, uint64_t sigma,
        std::vector<uint64_t>& currcode) {}


    tunstall_coder() {
      D_size = 65536;
    }

    tunstall_coder(std::vector<uint32_t> &seq, uint32_t block_size, uint64_t D_size_init) {

      uint64_t i;
      D_size = 1 << w;

      bSize = block_size;

      //std::cout << "Create alphabet..." << endl;
      map<uint32_t, uint32_t> alphabet_map;

      for (i = 0; i < seq.size(); ++i)
        alphabet_map[seq[i]] = 1;

      sigma = alphabet_map.size();

      i = 0;
      for (map<uint32_t, uint32_t>::iterator it = alphabet_map.begin(); it != alphabet_map.end(); ++it) {
        it->second = i;
        i++;
        map_table.push_back(it->first);
      }

      //std::cout << "Freq Table..." << endl;
      vector<uint64_t> freq_table(sigma, 0);

      for (i = 0; i < seq.size(); ++i) {
        ++freq_table[alphabet_map[seq[i]]];
      }

      uint64_t Pmin = freq_table[0];

      for (i = 1; i < freq_table.size(); ++i)
        if (freq_table[i] < Pmin)
          Pmin = freq_table[i];

      //std::cout << "Heap..." << endl;
      priority_queue<heap_node> H;

      vector<tree_node> tree(sigma);

      // tree initially contains sigma elements,
      // and they have no children in the tree (indicated with -1)

      for (i = 0; i < sigma; ++i) {
        tree[i] = tree_node(0, -1);
        H.push(heap_node((float)freq_table[i] / seq.size(), i));
      }

      // now, tree nodes are expanded according to their probabilities

      // std::cout << "Expansion..." << endl;
      float probMinLeaf = (float)Pmin / freq_table.size();
      uint dictionary_size = sigma;
      while (dictionary_size + sigma <= D_size) {
        heap_node N = H.top();
        H.pop();
        tree[N.second].second = tree.size();  // pointer to the children

        for (i = 0; i < sigma; ++i) {
          pair<float, uint64_t> p(N.first * ((float)freq_table[i] / seq.size()), tree.size());
          H.push(p);

          tree.push_back(tree_node(0, -1));
        }
        dictionary_size = dictionary_size + sigma;
      }

      // Now, traverse the tree to store the codes into the dictionary
      uint64_t curindex = 0;

      //std::cout << "Traverse..." << endl;
      D = std::vector<vector<uint64_t>>(D_size);

      vector<uint64_t> currcode;

      traverse(tree, 0, curindex, sigma, currcode);

      //std::cout << "Prefix sum..." << endl;
      uint64_t curnode = 0;
      uint64_t prefix_sum = 0, nelems_block = 0;

      blockElement bElem;

      bElem.prefix_sum = 0;
      bElem.starting_position = 0;

      block.push_back(bElem);

      std::vector<uint32_t> compressed_seq_aux;
      for (i = 0; i < seq.size(); ++i) {
        prefix_sum += seq[i];
        nelems_block++;
        if (tree[curnode + alphabet_map[seq[i]]].second != -1) {
          if (nelems_block == block_size) {
            compressed_seq_aux.push_back(tree[curnode + alphabet_map[seq[i]]].first);
            bElem.prefix_sum = prefix_sum;
            bElem.starting_position = compressed_seq_aux.size();
            block.push_back(bElem);
            nelems_block = 0;
            curnode = 0;
            nelems_block = 0;
          } else
            curnode = tree[curnode + alphabet_map[seq[i]]].second;
        } else {
          compressed_seq_aux.push_back(tree[curnode + alphabet_map[seq[i]]].first);
          curnode = 0;  // go back to the Tunstall tree root again
          if (nelems_block == block_size) {
            bElem.prefix_sum = prefix_sum;
            bElem.starting_position = compressed_seq_aux.size();
            block.push_back(bElem);
            nelems_block = 0;
          }
        }
      }

      if (nelems_block > 0)
        compressed_seq_aux.push_back(tree[curnode + alphabet_map[seq[seq.size() - 1]]].first);

      //std::cout << "Copying to a compressed_seq int_vector..." << "\n";
      compressed_seq = sdsl::int_vector<w>(compressed_seq_aux.size());
      for(uint64_t i = 0; i < compressed_seq_aux.size(); i++) {
        compressed_seq[i] = compressed_seq_aux[i];
      }
    }

    uint64_t decode(uint64_t i) {
      uint64_t b_ = i / bSize;

      uint64_t sum = block[b_].prefix_sum;
      uint64_t p = block[b_].starting_position;

      uint64_t j, size, nDecode = i % bSize + 1;

      for (j = 0; j <= nDecode; ++p) {
        size = D[compressed_seq[p]].size();
        if (j + size < nDecode) {
          sum += D[compressed_seq[p]][size - 1];
          j += size;
        } else {
          sum += D[compressed_seq[p]][nDecode - j - 1];
          break;
        }
      }

      return sum;
    }
    uint64_t dict_size() {
      // Tunstall dictionary size, in bytes  
      uint64_t i;
      uint64_t total_size = 0;
      for (i = 0; i < D.size(); ++i) {
        total_size += sizeof(uint16_t) * D[i].size();
      }
      return total_size;
    }
    uint64_t compressed_seq_size() { 

      return compressed_seq.bit_size() / 8;
      //return sdsl::size_in_bytes(compressed_seq);
    }

    uint64_t block_vec_size() {
      return block.size() * sizeof(blockElement);
    }

    uint64_t size() {
      // compressed size, in bytes 
      return dict_size() + compressed_seq_size() + block_vec_size();
    }
    uint64_t nCodewords() {
      return compressed_seq.size();
    }
};

#endif
