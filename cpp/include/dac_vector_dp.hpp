/*! \file dac_vector.hpp
   \brief dac_vector.hpp contains a vector which stores the values with variable length codes.
   \author Simon Gog, Niklas Baumstark, Julian Labeit
*/
#ifndef SDSL_DAC_VECTOR
#define SDSL_DAC_VECTOR

//#include <sdsl/bits.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/iterators.hpp>
#include <sdsl/rank_support_v5.hpp>
#include <sdsl/rrr_vector.hpp>

//! Namespace for the succinct data structure library.
namespace sdsl
{

//! A generic immutable space-saving vector class for unsigned integers.
/*! The values of a dac_vector are immutable after the constructor call.
 *  The ,,escaping'' technique is used to encode values. Bit widths of each
 *  encoding level are chosen optimally via dynamic programming.
 *  \par References
 *       [1] N. Brisaboa and S. Ladra and G. Navarro: ,,DACs: Bringing Direct Access to Variable-Length Codes'',
             Information Processing and Management (IPM) 2013
 *
 * \tparam t_bv                    Bit vector to use for overflow bits. Use
 *                                 rrr_vector<> for maximum compression, and
 *                                 bit_vector for speed.
 * \tparam t_default_max_levels    Maximum number of levels to use.
 */
template <typename t_bv = bit_vector, int t_default_max_levels = 64>
class dac_vector_dp
{
        static_assert(t_default_max_levels > 0, "invalid max level count");
    public:
        typedef typename int_vector<>::value_type        value_type;
        typedef random_access_const_iterator<dac_vector_dp> const_iterator;
        typedef t_bv                                     overflow_bv;
        typedef typename t_bv::rank_1_type               overflow_bv_rank1;
        typedef const_iterator                           iterator;
        typedef const value_type                         const_reference;
        typedef const_reference                          reference;
        typedef const_reference*                         pointer;
        typedef const pointer                            const_pointer;
        typedef int_vector<>::size_type                  size_type;
        typedef ptrdiff_t                                difference_type;
    private:
        size_t m_size;
        bit_vector m_overflow_tmp;
        overflow_bv m_overflow;
        overflow_bv_rank1 m_overflow_rank;
        std::vector<int_vector<>> m_data;
        std::vector<size_t> m_offsets;

        template <typename Container>
        void construct_level(
                size_t level, size_t overflow_offset,
                const std::vector<int>& bit_sizes,
                Container&& c)
        {
            if (level == bit_sizes.size()) {
                assert(c.size() == 0);
                assert(overflow_offset == m_overflow_tmp.size());
                // pack overflow bit vector (no-op if t_bv = bit_vector)
                m_overflow = overflow_bv(std::move(m_overflow_tmp));
                m_overflow_rank = overflow_bv_rank1(&m_overflow);
                m_overflow_tmp = bit_vector();
                return;
            }

            m_offsets.push_back(overflow_offset);
            size_t n = c.size();

            // mark elements with < *bit_sizes bits
            int bits_next = bit_sizes[level];
            size_t overflows = 0;
            int max_msb = 0; // max MSB of all values in c

            for (size_t i = 0; i < n; ++i) {
                int msb = bits::hi(c[i]);
                max_msb = std::max(max_msb, msb);
                if (msb >= bits_next) {
                    m_overflow_tmp[overflow_offset + i] = 1;
                    overflows++;
                } else {
                    m_overflow_tmp[overflow_offset + i] = 0;
                }
            }

            auto& data = m_data[level];
            data = int_vector<>(n - overflows);

            size_t idx_data = 0, idx_recurse = 0;
            int_vector<> recurse(overflows, 0, max_msb + 1);
            for (size_t i = 0; i < n; ++i) {
                if (m_overflow_tmp[overflow_offset + i]) {
                    recurse[idx_recurse++] = c[i];
                } else {
                    data[idx_data++] = c[i];
                }
            }
            util::bit_compress(data);
            assert(idx_data == n-overflows);
            assert(idx_recurse == overflows);

            construct_level(
                level + 1,
                overflow_offset + n,
                bit_sizes,
                recurse);
        }

    public:
        // copy-and-swap
        dac_vector_dp() = default;
        dac_vector_dp(const dac_vector_dp& other)
            : m_size(other.m_size)
            , m_overflow(other.m_overflow)
            , m_overflow_rank(other.m_overflow_rank)
            , m_data(other.m_data)
            , m_offsets(other.m_offsets)
        {
            m_overflow_rank.set_vector(&m_overflow);
        }

        void swap(dac_vector_dp& other) {
            std::swap(m_size, other.m_size);
            m_overflow.swap(other.m_overflow);
            util::swap_support(m_overflow_rank, other.m_overflow_rank,
                               &m_overflow, &other.m_overflow);
            std::swap(m_data, other.m_data);
            std::swap(m_offsets, other.m_offsets);
        }

        dac_vector_dp(dac_vector_dp&& other) : dac_vector_dp() {
            this->swap(other);
        }

        dac_vector_dp& operator=(dac_vector_dp other) {
            this->swap(other);
            return *this;
        }

        double cost(size_t n, size_t m) {
            // FIXME(niklasb) this seems to work well for both plain and rrr, but
            // can probably be improved
            double overhead = 128;
            if (n == 0 || m == 0 || m == n) return overhead;
            double plain = 1.02 * n;
            double entropy = (1.*m/n * log(1.*n/m) / log(2) +
                1.*(n-m)/n * log(1.*n/(n-m)) / log(2));
            double rrr = overhead + (0.1 + entropy) * n;
            return std::min(plain, rrr);
        }

        //! Constructor for a Container of unsigned integers.
        /*! \param c          A container of unsigned integers.
          * \param max_level  Maximum number of levels to use.
          */
        template<class Container>
        dac_vector_dp(Container&& c, int max_levels = t_default_max_levels) {
            assert(max_levels > 0);
            m_size = c.size();
            std::vector<uint64_t> cnt(128, 0);
            cnt[0] = m_size;
            int max_msb = 0;
            for (size_t i = 0; i < m_size; ++i) {
                auto x = c[i] >> 1;
                int lvl = 1;
                while (x > 0) {
                    cnt[lvl] += 1;
                    max_msb = std::max(max_msb, lvl);
                    x >>= 1;
                    ++lvl;
                }
            }

            // f[i][j] = minimum cost for subsequence with MSB >= i, when we can
            // use up to j levels.
            double f[max_msb + 2][max_levels + 1];
            int nxt[max_msb + 2][max_levels + 1];
            std::fill(f[max_msb + 1], f[max_msb + 1] + max_levels + 1, 0.0);
            std::fill(nxt[max_msb + 1], nxt[max_msb + 1] + max_levels + 1, -1);
            for (int b = max_msb; b >= 0; --b) {
                std::fill(f[b], f[b] + max_levels + 1,
                    std::numeric_limits<double>::infinity());
                for (int lvl = 1; lvl <= max_levels; ++lvl) {
                    for (int b2 = b+1; b2 <= max_msb + 1; ++b2) {
                        double w = b2*(cnt[b] - cnt[b2]) + cost(cnt[b], cnt[b2]) + f[b2][lvl - 1];
                        if (w < f[b][lvl]) {
                            f[b][lvl] = w;
                            nxt[b][lvl] = b2;
                        }
                    }
                }
            }
            std::vector<int> bit_sizes;
            int b = 0, lvl = max_levels;
            while (nxt[b][lvl] != -1) {
                b = nxt[b][lvl];
                lvl--;
                bit_sizes.push_back(b);
            }
            assert(bit_sizes.size() <= max_levels);

            size_t total_overflow_size = 0;
            for (size_t i = 0; i < c.size(); ++i) {
                size_t b = 0;
                int msb = bits::hi(c[i]);
                ++total_overflow_size;
                while (b < bit_sizes.size() && msb >= bit_sizes[b]) {
                    ++b;
                    ++total_overflow_size;
                }
            }

            m_data.resize(bit_sizes.size());
            m_overflow_tmp.resize(total_overflow_size);
            construct_level(0, 0, bit_sizes, c);
        }

        //! The number of levels in the dac_vector.
        size_t levels() const {
            return m_data.size();
        }

        //! The number of elements in the dac_vector.
        size_type size() const {
            return m_size;
        }

        //! Returns if the dac_vector is empty.
        bool empty() const { return !size(); }

        //! Iterator that points to the first element of the dac_vector.
        const const_iterator begin() const
        {
            return const_iterator(this, 0);
        }


        //! Iterator that points to the position after the last element of the dac_vector.
        const const_iterator end()const
        {
            return const_iterator(this, size());
        }

        //! []-operator
        value_type operator[](size_type i)const
        {
            size_t level = 0, offset = m_offsets[level];
            while (m_overflow[offset + i]) {
                i = m_overflow_rank(offset + i) - m_overflow_rank(offset);
                level++;
                offset = m_offsets[level];
            }
            i -= m_overflow_rank(offset + i) - m_overflow_rank(offset);
            return m_data[level][i];
        }

        //! Serializes the dac_vector to a stream.
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr,
                std::string name="") const {
            structure_tree_node* child = structure_tree::add_child(
                                            v, name, util::class_name(*this));
            size_type written_bytes = 0;

            written_bytes += m_overflow.serialize(out, child, "overflow");
            written_bytes += m_overflow_rank.serialize(out, child, "overflow_rank");

            written_bytes += sdsl::serialize(m_data, out, child, "data");
            written_bytes += sdsl::serialize(m_offsets, out, child, "offsets");

            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Load from a stream.
        void load(std::istream& in) {
            m_overflow.load(in);
            m_overflow_rank.load(in);
            m_overflow_rank.set_vector(&m_overflow);
            sdsl::load(m_data, in);
            sdsl::load(m_offsets, in);
            m_size = 0;
            for (auto& v : m_data)
                m_size += v.size();
        }
};
}
#endif
