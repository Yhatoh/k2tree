#include "huffman_coder.hpp"

using namespace std;

uint32_t MAX(uint32_t a, uint32_t b) {
    return  a < b ? b : a;
}

huffman_coder::huffman_coder() {;}

huffman_coder::~huffman_coder() {
    delete [] syms;
}


uint64_t huffman_coder::size_bits() {
    uint64_t size_lut =  LUT_SIZE * sizeof(uint32_t*);
    uint64_t size_compressed_seq = compressed_seq.size()*sizeof(uint32_t);
    uint64_t size_block_info = sdsl::size_in_bytes(block_info_starting_position) + sdsl::size_in_bytes(block_info_prefix_sum);
    uint64_t size_syms = n_symbols*ceil(log2(n_symbols)); // in bits
    uint64_t size_small_tables = 4*L*sizeof(uint32_t);

    return 8*(size_lut+size_compressed_seq+size_block_info+size_small_tables) + size_syms;
}

float huffman_coder::size_bpe()
{
    uint64_t size_lut =  LUT_SIZE * sizeof(uint32_t*);
    uint64_t size_compressed_seq = compressed_seq.size()*sizeof(uint32_t);
    uint64_t size_block_info = sdsl::size_in_bytes(block_info_starting_position) + sdsl::size_in_bytes(block_info_prefix_sum);
    uint64_t size_syms = n_symbols*ceil(log2(n_symbols)); // in bits
    uint64_t size_small_tables = 4*L*sizeof(uint32_t);
    
    return ((float)(8*(size_lut+size_compressed_seq+size_block_info+size_small_tables) + size_syms))/length;
}

uint64_t huffman_coder::size()
{
    uint64_t size_lut = LUT_SIZE * sizeof(uint32_t*);
    uint64_t size_compressed_seq = compressed_seq.size()*sizeof(uint32_t);
    uint64_t size_block_info = sdsl::size_in_bytes(block_info_starting_position) + sdsl::size_in_bytes(block_info_prefix_sum);
    uint64_t size_syms = n_symbols*ceil(log2(n_symbols))/8;
    uint64_t size_small_tables = 4*L*sizeof(uint32_t);
    
    return size_lut+size_compressed_seq+size_block_info+size_small_tables + size_syms;
}

uint32_t huffman_coder::sigma() {return n_symbols;};

uint64_t huffman_coder::nCodewords() { return compressed_seq.size(); }

void huffman_coder::build_freq_table(vector<uint32_t> &seq)
{
    uint64_t i;
    
    freq_table = new uint32_t[max_symbol + 2];
    freq_table[0] = 1;  // EOF, not needed in the next version
    
    for (i = 1; i < max_symbol+2; ++i)
        freq_table[i] = 0;

    for (i = 0; i < seq.size(); ++i)
        freq_table[seq[i]]++;
}


/*
** Build lj_base[] and offset from the codelens in A[0..n-1]
** A[] need not be sorted.
**
** Return cw_lens[] a freq count of codeword lengths.
*/
void huffman_coder::build_canonical_arrays(uint32_t max_cw_length)
{
    uint32_t* q;
    uint32_t* p;

    // build offset
    q = offset;
    *q = 0;
    for (p = cw_lens + 1; p < cw_lens + max_cw_length; p++, q++)
        *(q + 1) = *q + *p;

    // generate the min_code array
    q = min_code + max_cw_length - 1;
    *q = 0;
    for (q--, p = cw_lens + max_cw_length; q >= min_code; q--, p--)
        *q = (*(q + 1) + *p) >> 1;

    // generate the lj_base array
    q = lj_base;
    uint32_t* pp = min_code;
    int32_t left_shift = (sizeof(uint32_t) << 3) - 1;
    for (p = cw_lens + 1; q < lj_base + max_cw_length;
         p++, q++, pp++, left_shift--)
        if (*p == 0)
            *q = *(q - 1);
        else
            *q = (*pp) << left_shift;
    for (p = cw_lens + 1, q = lj_base; *p == 0; p++, q++)
        *q = MAX_ULONG;

} // build_canonical_arrays()


/*
** INPUT: syms[0..n-1] lists symbol numbers
**        freq[i] contains the codeword length of symbol i
**        cw_lens[1..max_cw_length] is the number of codewords of length i
**
** OUTPUT: None
**
** SIDE EFFECTS: syms[0..max_symbol] is overwritten with canonical code mapping.
**               cw_lens[] is destroyed.
*/
void huffman_coder::generate_mapping(uint32_t max_cw_length, uint32_t n)
{
    int32_t i;

    for (i = 1; i <= (int)max_cw_length; i++)
        cw_lens[i] += cw_lens[i - 1];

    for (i = n - 1; i >= 0; i--) {
    	  lens[i] = syms[syms[i]] = cw_lens[freq_table[syms[i]] - 1]++;    
    }

} /* generate_mapping() */


void huffman_coder::calculate_minimum_redundancy(int32_t n)
{
    int32_t root; /* next root node to be used */
    int32_t leaf; /* next leaf to be used */
    int32_t next; /* next value to be assigned */
    int32_t avbl; /* number of available nodes */
    int32_t used; /* number of internal nodes */
    uint32_t dpth; /* current depth of leaves */

    /* check for pathological cases */
    if (n == 0) {
        return;
    }
    if (n == 1) {
        freq_table[syms[0]] = 0;
        return;
    }

    /* first pass, left to right, setting parent pointers */
    freq_table[syms[0]] += freq_table[syms[1]];
    root = 0;
    leaf = 2;
    for (next = 1; next < n - 1; next++) {
        /* select first item for a pairing */
        if (leaf >= n || freq_table[syms[root]] < freq_table[syms[leaf]]) {
            freq_table[syms[next]] = freq_table[syms[root]];
            freq_table[syms[root++]] = next;
        } else
            freq_table[syms[next]] = freq_table[syms[leaf++]];

        /* add on the second item */
        if (leaf >= n || (root < next && freq_table[syms[root]] < freq_table[syms[leaf]])) {
            freq_table[syms[next]] += freq_table[syms[root]];
            freq_table[syms[root++]] = next;
        } else
            freq_table[syms[next]] += freq_table[syms[leaf++]];
    }

    /* second pass, right to left, setting internal depths */
    freq_table[syms[n - 2]] = 0;
    for (next = n - 3; next >= 0; next--)
        freq_table[syms[next]] = freq_table[syms[freq_table[syms[next]]]] + 1;

    /* third pass, right to left, setting leaf depths */
    avbl = 1;
    used = dpth = 0;
    root = n - 2;
    next = n - 1;
    while (avbl > 0) {
        while (root >= 0 && freq_table[syms[root]] == dpth) {
            used++;
            root--;
        }
        while (avbl > used) {
            freq_table[syms[next--]] = dpth;
            avbl--;
        }
        avbl = 2 * used;
        dpth++;
        used = 0;
    }
}


int32_t pcmp(char* a, char* b) { return *((uint32_t*)a) - *((uint32_t*)b); }

void huffman_coder::build_codes()
{
    uint32_t i, *p;
    uint32_t max_codeword_length; //, min_codeword_length;

    indirect_sort(freq_table, syms, syms, n_symbols);

    calculate_minimum_redundancy(n_symbols);

    // calculcate max_codeword_length and set cw_lens[]
    for (i = 0; i <= L; i++)
        cw_lens[i] = 0;
    // min_codeword_length = max_codeword_length = freq[syms[0]];
    max_codeword_length = 0;
    for (p = syms; p < syms + n_symbols; p++) {
        if (freq_table[*p] > max_codeword_length)
            max_codeword_length = freq_table[*p];
        cw_lens[freq_table[*p]]++;
    }

    build_canonical_arrays(max_codeword_length);

    nqsort((char*)syms, n_symbols, sizeof(uint32_t), pcmp);

    generate_mapping(max_codeword_length, n_symbols);
    
    max_cw_len = max_codeword_length;

} /* build_codes() */


void huffman_coder::build_lut()
{
    uint32_t max, min; // range of left justified "i"
    int32_t i, j = max_cw_len - 1; // pointer into lj

    for (i = 0; i < LUT_SIZE; i++) {
        min = i << ((sizeof(uint32_t) << 3) - LUT_BITS);
        max = min | MAX_IT;

        while ((j >= 0) && (max > lj_base[j]))
            j--;

        // we know max is in range of lj[j], so check min
        if (min >= lj_base[j + 1])
            lut[i] = lj_base + j + 1;
        else
            lut[i] = NULL; //-(j+1);
    }

} /* build_lut() */


void huffman_coder::OUTPUT_ULONG(uint32_t n, char len)
{
    if ((uint32_t)len < buff_btg) {
        buff <<= len;
        buff |= n;
        buff_btg -= len;
    } else {
        buff <<= buff_btg;
        buff |= (n) >> (len - buff_btg);
        compressed_seq.push_back(buff);
        buff = n & ~((~0) << (len - buff_btg)); // OJO, verificar esto!!!
        buff_btg = BUFF_BITS - (len - buff_btg);
    }
}


//
// Interpret the next len bits of the input as a ULONG and return the result
//
uint32_t huffman_coder::INPUT_ULONG(uint32_t &cur_int, int32_t len)
{
    if (len == 0)
        return 0;

    uint32_t n;

    if (buff_btg == BUFF_BITS)
        n = (buff) >> (BUFF_BITS - len);
    else
        n = ((buff) << (BUFF_BITS - buff_btg)) >> (BUFF_BITS - len);

    if ((uint32_t)len < buff_btg)
        buff_btg -= len;
    else {
        len -= buff_btg;
        cur_int++;
        buff = compressed_seq[cur_int];        
        buff_btg = BUFF_BITS;
        if (len > 0) {
            n |= (buff) >> (BUFF_BITS - len);
            buff_btg -= len;
        }
    }

    if (buff_btg == 0) {
        buff = compressed_seq[cur_int];
        cur_int++;
        buff_btg = BUFF_BITS;
    }
    
    return n;
} // INPUT_ULONG()


/*
** Canonical encode.  cwlens[] contains codeword lens, mapping[] contains
** ordinal symbol mapping.
*/
uint32_t huffman_coder::output(uint32_t i, uint32_t mapping[], uint32_t cwlens[])
{
    uint32_t sym_num = mapping[i]; // ordinal symbol number

    uint32_t len = cwlens[i];

    uint32_t cw = min_code[len - 1] + (sym_num - offset[len - 1]);

    OUTPUT_ULONG(cw, len);

    return len;
} /* output() */


void huffman_coder::encode(vector<uint32_t> &seq, uint32_t _block_size) 
{    
    uint64_t n, total_len;
    uint64_t i;

    length = seq.size();
            
    map<uint32_t, uint32_t> alphabet_map;
    max_symbol = seq[0]+1;    
    for (i = 0; i < seq.size(); i++) {
        alphabet_map[seq[i]] = 1;
        seq[i] = seq[i] + 1;        
        if (seq[i] > max_symbol)
            max_symbol = seq[i];
    }

    n_symbols = alphabet_map.size()+1;

    alphabet_map.clear();
    
    lens = new int32_t[n_symbols+1];    
      
    syms = new uint32_t[max_symbol + 2];
    
    uint32_t* mapping = new uint32_t[max_symbol + 2];
    
    build_freq_table(seq);

    freq_table[0] = 1;
    n = 0;
    for (i = 0; i <= max_symbol; i++)
        if (freq_table[i] > 0)
            syms[n++] = i;
            
    for (i = 0; i < max_symbol; i++)
        mapping[i] = syms[i];    
    
    build_codes();
    
    for (min_cw_len = 0; cw_lens[min_cw_len] == 0; min_cw_len++)
        ;
    min_cw_len++; // CHECK THIS!

    block_size = _block_size;

    buff = 0;
    buff_btg = BUFF_BITS; 
    
    uint64_t sum = 0;

    sdsl::int_vector<> block_pos(seq.size()/block_size + 1), block_sum(seq.size()/block_size + 1);
     
    uint64_t last_block = 0;
    
    for (i = 0; i < seq.size(); ++i) {
        if (i%block_size == 0) {
            if (i != 0) compressed_seq.push_back(buff<<buff_btg);
            buff = 0;
            buff_btg = BUFF_BITS;        
            
            block_pos[last_block] = compressed_seq.size();            
            block_sum[last_block++] = sum;
            
        }
        sum += (seq[i] - 1);  // -1 because of the +1 done to the original sequence
        total_len += output(seq[i], syms, freq_table);
    }
    
    block_info_starting_position = sdsl::enc_vector<>(block_pos);
    block_info_prefix_sum = sdsl::enc_vector<>(block_sum);
    
    compressed_seq.push_back(buff<<buff_btg);
    
    build_lut();
    
    
    n = n_symbols;

    for (i = 0; i <= max_symbol; i++)
        syms[i] = mapping[i];

    delete [] mapping;
    
    uint32_t t, from, S;
    int32_t start = 0;    
    lens[n] = 1; // sentinel
    while (start < (int32_t)n) {
        from = start;
        S = syms[start];

        while (lens[from] >= 0) {
            i = lens[from];
            lens[from] = -1;
            t = syms[i];
            syms[i] = S;
            S = t;
            from = i;
        }

        while (lens[start] == -1)
            start++; // find next start (if any)
    }

    delete lens;
    delete freq_table;
    
}


uint64_t huffman_coder::decode(uint64_t i)
{
    uint32_t n_decoded, cur_int;
    uint32_t code = 0;
    uint32_t bits_needed = sizeof(uint32_t) << 3;
    uint32_t currcode;
    uint32_t currlen = sizeof(uint32_t) << 3;
    uint32_t* lj;
    uint32_t* start_linear_search = lj_base + MAX(LUT_BITS, min_cw_len) - 1;

    uint64_t sum = block_info_prefix_sum[i/block_size];
    cur_int = block_info_starting_position[i/block_size];
    i = i - (i/block_size)*block_size;

    buff = compressed_seq[cur_int];
    buff_btg = BUFF_BITS;

    for (n_decoded = 0; n_decoded <= i; ++n_decoded) {
        code |= INPUT_ULONG(cur_int, bits_needed);

        lj = lut[code >> ((sizeof(uint32_t) << 3) - LUT_BITS)];
        if (lj == NULL)
            for (lj = start_linear_search; code < *lj; lj++)
                ;
        currlen = lj - lj_base + 1;

        // calculate symbol number
        currcode = code >> ((sizeof(uint32_t) << 3) - currlen);
        currcode -= min_code[currlen - 1];
        currcode += offset[currlen - 1];

        sum += syms[currcode]-1;    
                    
        code <<= currlen;
        bits_needed = currlen;
    }
    
    return sum;
}  
