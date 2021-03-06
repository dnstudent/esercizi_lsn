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

/**
 * Converts from base 2^12 to base 10
 * @param a,b,c,d Base-2^12 digits of the number to be converted from the most to the least significant.
 */
template<typename uint>
constexpr uint b4096tob10(uint a, uint b, uint c, uint d) noexcept {
    return d + static_cast<uint>(4096U) *
                       (c + static_cast<uint>(4096U) * (b + static_cast<uint>(4096U) * a));
}

/**
 * Converts from base 10 to base 2^12.
 * @param x The decimal number to be converted. Must be < 2^48.
 * @param a, b, c, d Pointers to the base-2^12 digits from the most to the least significant.
 */
template<typename uint>
constexpr void b10tob4096(uint x, uint *a, uint *b, uint *c, uint *d) {
    if (x < 4096UL * 4096UL * 4096UL * 4096UL) { throw; }
    *d = x % 4096U;
    x /= 4096U;
    *c = x % 4096U;
    x /= 4096U;
    *b = x % 4096U;
    x /= 4096U;
    *a = x % 4096U;
}

/**
 * Class implementing the RANNYU routine, a 48-bit linear congruential generator with specific choices for the multiplicator and the prime addends.
 * Computation is performed in base 2^12.
 * The procedure was adapted from the provided file to fit with modern c++ algorithms and to be built in a cmake project.
 * Global variables were replaced with a class managing the procedure state and URBG requirements were implemented.
 */
class ARandom {
protected:
public:
    typedef uint_fast64_t result_type;
    // constructors
    ARandom() = default;
    explicit ARandom(result_type s) : ARandom() { seed(s); }
    ARandom(const std::string_view &, const std::string_view &, size_t = 1);

    /*------------------
     * URBG methods: they enable ARandom to be used with the stanadrd library <random> distributions
     ------------------*/
    result_type operator()();
    [[nodiscard]] [[maybe_unused]] static constexpr result_type min() { return 0UL; }
    [[nodiscard]] [[maybe_unused]] static constexpr result_type max() {
        // 2^48 - 1
        return 281474976710655UL;
    }

    /**
     * Seeds the generator with a base 10 number
     * @param s Seed in base 10
     */
    [[maybe_unused]] void seed(result_type s);

    /**
     * Sets the seed and the two least significant base-2^12 digits of the prime (unsigned decimals < 4096)
     * @param seed Seed in base 2^12 as four decimals
     * @param prime3 Second least significant base-2^12 digit
     * @param prime4 Least significant base-2^12 digit
     */
    void SetRandom(result_type const *seed, result_type prime3, result_type prime4);
    void SaveSeed() const;
    [[maybe_unused]] double Rannyu();
    [[maybe_unused]] double Rannyu(double min, double max);
    [[maybe_unused]] double Gauss(double mean, double sigma);

private:
    // multiplyer in base 2^12
    const result_type m_m1{502UL}, m_m2{1521UL}, m_m3{4071UL}, m_m4{2107UL};
    // most significant digits of the prime in base 2^12
    const result_type m_p1{0UL}, m_p2{0UL};
    // seed
    result_type m_l1{0UL}, m_l2{0UL}, m_l3{0UL}, m_l4{1UL};
    // least significant digits of the prime
    result_type m_p3{2892UL}, m_p4{2587UL};
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
