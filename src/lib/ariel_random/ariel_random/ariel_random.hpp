/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef ARIEL_RANDOM_
#define ARIEL_RANDOM_

#include <cstddef>
#include <limits>
#include <string_view>

template<typename uint>
constexpr uint b4096tob10(uint a, uint b, uint c, uint d) {
    return d + static_cast<uint>(4096U) *
                       (c + static_cast<uint>(4096U) * (b + static_cast<uint>(4096U) * a));
}

template<typename uint>
constexpr void b10tob4096(uint x, uint *a, uint *b, uint *c, uint *d) {
    *d = x % 4096U;
    x /= 4096U;
    *c = x % 4096U;
    x /= 4096U;
    *b = x % 4096U;
    x /= 4096U;
    *a = x % 4096U;
}

// TODO: cleanup
class ARandom {
protected:
public:
    typedef uint64_t result_type;
    // constructors
    ARandom() = default;
    ARandom(ARandom &x) = default;
    [[maybe_unused]] explicit ARandom(result_type s) : ARandom() { seed(s); }
    template<typename Sseq>
    [[maybe_unused]] explicit ARandom(Sseq &q) : ARandom() {
        seed(q);
    }
    ARandom(const std::string_view &, const std::string_view &, size_t = 1);
    [[maybe_unused]] void seed();
    void seed(result_type s);
    template<typename Sseq>
    void seed(Sseq & /*q*/) {}
    [[maybe_unused]] void discard(unsigned long long z) {
        for (auto i = 0ULL; i < z; i++) operator()();
    }
    inline bool operator==(const ARandom &oth) const {
        return (m_l1 == oth.m_l1 && m_l2 == oth.m_l2 && m_l3 == oth.m_l3 && m_l4 == oth.m_l4 &&
                m_m1 == oth.m_m1 && m_m2 == oth.m_m2 && m_m3 == oth.m_m3 && m_m4 == oth.m_m4 &&
                m_n1 == oth.m_n1 && m_n2 == oth.m_n2 && m_n3 == oth.m_n3 && m_n4 == oth.m_n4);
    }
    inline bool operator!=(const ARandom &oth) const { return !(*this == oth); }
    // Should be implemented but who cares
    template<typename os>
    friend os &operator<<(os & /*o*/, const ARandom & /*rng*/) {}
    template<typename is>
    friend is &operator>>(is & /*i*/, ARandom & /*rng*/) {}
    // methods
    [[maybe_unused]] void SetRandom(result_type const *, result_type, result_type);
    [[maybe_unused]] void SaveSeed() const;
    [[maybe_unused]] double Rannyu();
    result_type operator()();
    [[nodiscard]] [[maybe_unused]] static constexpr result_type min() { return 0UL; }
    [[nodiscard]] [[maybe_unused]] static constexpr result_type max() {
        // 2^48 - 1
        return 281474976710655UL;
    }
    [[maybe_unused]] double Rannyu(double min, double max);
    [[maybe_unused]] double Gauss(double mean, double sigma);

private:
    // moltiplicatore
    const result_type m_m1{502UL}, m_m2{1521UL}, m_m3{4071UL}, m_m4{2107UL};
    // prime due word del prime
    result_type m_n1{0UL}, m_n2{0UL};
    // seed
    result_type m_l1{0UL}, m_l2{0UL}, m_l3{0UL}, m_l4{1UL};
    // seconde due word del prime
    result_type m_n3{2892UL}, m_n4{2587UL};
};

#endif

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
