#ifndef __BIT_VECTOR__
#define __BIT_VECTOR__

// c++ includes
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "debug.hpp"
// super simple append only bit_vector
struct bvector {
  std::vector< uint8_t > bv;
  uint64_t n;

  bvector() {
    n = 0;
  }

  uint64_t size() {
    return n;
  }

  uint8_t operator[](const uint64_t i) {
    assert(i <= n);
    const uint64_t integer = i / 8;
    const uint64_t i_int = i % 8;

    return ((bv[integer] & (1 << i_int)) != 0);
  }

  void set(const uint64_t num, const uint64_t i, const uint64_t w) {
    assert(i + w <= n);
    uint64_t integer = i / 8;
    uint64_t i_int = i % 8;

    for(uint64_t bit = 0; bit < w; bit++) {
      if(num & (1 << bit)) {
        bv[integer] |= 1 << i_int;
      } else {
        bv[integer] &= ~(1 << i_int);
      }
      i_int++;
      if(i_int == 8) {
        integer++;
        i_int = 0;
      }
    }
  }
    
  void reserve(uint64_t i) {
    bv.reserve((i + 8 - 1) / 8);
  }

  void push_back(uint8_t bit) {
    if(n % 8 == 0) {//full
      bv.push_back(0);
    }
    const uint64_t integer = n / 8;
    const uint64_t i_int = n % 8;
    n++;
    bv[integer] |= (bit << i_int);
  }

  void destroy() {
    std::vector< uint8_t >().swap(bv);
  }

  // this functions is super slow
  void concat(bvector& a, uint64_t i = 0, uint64_t end = 0) {
    assert(i <= a.size());
    assert(end <= a.size());
    if(end == 0) end = a.size();
    for(uint64_t j = i; j < end; j++) {
      push_back(a[j]);
    }
  }

  void clear() {
    bv.clear();
    n = 0;
  }
};

// a more complex append only bit vector with left to right scan

const size_t mod = (1 << 3) - 1;
const size_t size_int = 8;

struct block_vector {
  std::vector< std::vector< uint8_t > > blocks;
  std::vector< size_t > size_blocks;
  size_t curr_bit;
  // for left-right scan
  size_t pos, block, pos_block;
  size_t save_pos, save_block, save_pos_block;

  block_vector() : curr_bit(0), blocks({{}}), size_blocks({0}) {}

  inline size_t size() { return curr_bit; }
  inline void push_back(uint8_t bit) {
    curr_bit++;
    if((size_blocks.back() & mod) == 0) {
      blocks.back().push_back(0);
    }

    blocks.back()[size_blocks.back() / (mod + 1)] |= (bit << (size_blocks.back() & mod));
    size_blocks.back()++;
  }

  inline void concat(block_vector& bv) {
    blocks.reserve(blocks.size() + bv.blocks.size());
    blocks.insert(blocks.end(), bv.blocks.begin(), bv.blocks.end());

    size_blocks.reserve(size_blocks.size() + bv.size_blocks.size());
    size_blocks.insert(size_blocks.end(), bv.size_blocks.begin(), bv.size_blocks.end());

    curr_bit += bv.curr_bit;
  }

  inline void init_scan() {
    pos = block = pos_block = 0;
  }

  inline void move() {
    assert(pos < curr_bit);
    pos_block++;

    if(pos_block == size_blocks[block]) {
      pos_block = 0;
      block++;
    }
    pos++;
  }

  inline uint8_t read() {
    assert(pos < curr_bit);

    return (blocks[block][pos_block / (mod + 1)] & (1 << (pos_block & mod))) > 0;
  }

  inline uint8_t read(size_t block, size_t pos_block) {
    return (blocks[block][pos_block / (mod + 1)] & (1 << (pos_block & mod))) > 0;
  }

  inline void destroy() {
    std::vector< std::vector< uint8_t > >().swap(blocks);
    std::vector< size_t >().swap(size_blocks);
  }

  inline void swap(block_vector& bv) {
    std::swap(blocks, bv.blocks);
    std::swap(size_blocks, bv.size_blocks);
    std::swap(curr_bit, bv.curr_bit);
  }

  inline void save() {
    save_pos = pos;
    save_block = block;
    save_pos_block = pos_block;
  }

  // append from [save_pos, pos]
  inline void append_in(block_vector &bv) {
    if(block != save_block) {
      // copying several blocks
      for(; ((save_pos_block & mod) != 0) && (save_pos_block < size_blocks[save_block]); save_pos_block++) {
        bv.push_back(read(save_block, save_pos_block));
      }

      if(bv.blocks.back().size() == 0) {
        bv.blocks.pop_back();
        bv.size_blocks.pop_back();
      }

      if(save_pos_block < size_blocks[save_block]) {
        // now you are allow to just copy
        // from save_pos_block / (mod + 1)  to save_pos_block + size_blocks[save_block] / (mod + 1)
        bv.blocks.push_back({});
        bv.blocks.back().insert(bv.blocks.back().end(),
                                blocks[save_block].begin() + save_pos_block / (mod + 1),
                                blocks[save_block].begin() + (size_blocks[save_block]) / (mod + 1) + 1); // i hope

        bv.size_blocks.push_back(size_blocks[save_block] - save_pos_block);
      }
      save_block++;
      save_pos_block = 0;

      for(; save_block < block; save_block++) {
        bv.blocks.push_back({});
        bv.blocks.back().insert(bv.blocks.back().end(),
                                blocks[save_block].begin(),
                                blocks[save_block].end());
        bv.size_blocks.push_back(size_blocks[save_block]);
      }

      // last block to add
      bv.blocks.push_back({});
      bv.blocks.back().insert(bv.blocks.back().end(),
                              blocks[save_block].begin(),
                              blocks[save_block].begin() + (pos_block / (mod + 1)));
      bv.size_blocks.push_back((mod + 1) * (pos_block / (mod + 1)));
      save_pos_block += (mod + 1) * (pos_block / (mod + 1));
      for(; save_pos_block <= pos_block; save_pos_block++)
        bv.push_back(read(block, save_pos_block));

      bv.curr_bit += (pos - save_pos + 1);
      return;
    }

    for(; ((save_pos_block & mod) != 0) && (save_pos_block != pos_block); save_pos_block++) {
      bv.push_back(read(save_block, save_pos_block));
    }

    if(save_pos_block <= pos_block) {
      if(bv.blocks.back().size() == 0) {
        bv.blocks.pop_back();
        bv.size_blocks.pop_back();
      }
      bv.blocks.push_back({});
      bv.blocks.back().insert(bv.blocks.back().begin(),
                              blocks[block].begin() + save_pos_block / (mod + 1),
                              blocks[block].begin() + pos_block / (mod + 1));
      save_pos_block += (mod + 1) * (bv.blocks.back().size());
      bv.size_blocks.push_back((mod + 1) * (bv.blocks.back().size()));
      for(; save_pos_block <= pos_block; save_pos_block++)
        bv.push_back(read(block, save_pos_block));
    }

    bv.curr_bit += (pos - save_pos + 1);
    return;
  }
};



#endif // !__BIT_VECTOR__

