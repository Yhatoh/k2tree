#include <iostream>
#include <numeric>
#include <random>
#include <sdsl/bp_support_sada.hpp>
#include <sdsl/io.hpp>
#include <string>
#include <vector>
#include <stack>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/int_vector.hpp>

#include "bp_sdsl_idems.hpp"

template<typename K>
std::vector<K> generate_uniform(const size_t n, const K u, const int seed = 420) {
    
    static_assert(std::is_integral<K>::value);

    std::mt19937 gen(seed); 
    
    std::uniform_int_distribution<K> dis(1, u); 

    std::vector<K> results;
    results.reserve(n);

    for(auto i = 0; i < n; ++i) {
        results.push_back(dis(gen));
    }
    
    return results;
}

template<typename K>
sdsl::bit_vector build_cartesian_tree(const std::vector<K> &data) {
    
    sdsl::bit_vector::value_type max_excess = 0, curr_excess = 0;
    sdsl::bit_vector bp(data.size() * 2 + 2, 0);

    if(data.size() > 0) [[likely]] {
        size_t curr = 0, bp_curr = 0;
        std::stack<K> s;
        
        s.push(std::numeric_limits<K>::min());

        bp[bp_curr] = 1;
        bp_curr++;

        while(curr < data.size()) {

            while((data[curr] < s.top()) && s.size() > 1) {
                s.pop();
                bp_curr++;
                curr_excess--;
            }

            bp[bp_curr] = 1;
            curr_excess++;
            bp_curr++;

            if(curr_excess > max_excess) max_excess = curr_excess;

            s.push(data[curr]);

            curr++;
        }
    }

    // just for debugging
    std::cout << "tree height: " << max_excess << std::endl;

    return bp;  
}

int main() {

    const int32_t u = 1e9;
    const size_t n = 1e6; 

    std::vector<int32_t> data = generate_uniform<int32_t>(n, u);
    //std::vector<int32_t> data = {5, 2, 3, 1, 4, 3, 1, 4, 1, 2};

    sdsl::bit_vector bp = build_cartesian_tree<int32_t>(data);
    bp_support_sada<> x(&bp);
    cout << size_in_bytes(bp) * 8 + size_in_bytes(x) * 8 << endl;


    //for(const auto &b : bp)
    //    if(b) std::cout << "(";
    //    else std::cout << ")";
    std::cout << "----" << std::endl;

    bp_sdsl_idems<> attempt(bp);
    cout << attempt.size_in_bits() << " " << (double) attempt.size_in_bits() / attempt.nodes() << "\n";


    return 0;
}
