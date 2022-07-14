//
// Created by Davide Nicoli on 07/07/22.
//

#ifndef ESERCIZI_LSN_UNIFORM_INT_HPP
#define ESERCIZI_LSN_UNIFORM_INT_HPP

#include <cmath>
#include <random>
#include <type_traits>

#include <ariel_random/ariel_random.hpp>

namespace distributions {
    /**
     * An uniform int distribution to sample in [a, b]. For AppleClang this wraps floor(uniform_real_distribution) or floor(rng.Rannyu), as Apple libc++
     * apparently has a bug which prevents std::uniform_int_distribution procedure to work properly with ARandom.
     * Other OSes/compilers/standard libraries seem unaffected, so for them this is only a wrapper around std::uniform_int_distribution
     */
    template<typename IntType>
    class uniform_int {
        static_assert(std::is_integral_v<IntType>);

    public:
        typedef IntType result_type;
        uniform_int(IntType a, IntType b)
#ifdef __apple_build_version__
            : m_a{double(a)}, m_b{double(b + 1)}, m_d(m_a, m_b)
#else
            : m_id(a, b)
#endif
        {
        }

#ifdef __apple_build_version__
        inline result_type operator()(ARandom &rng) {
            return static_cast<result_type>(std::floor(rng.Rannyu(m_a, m_b)));
        }
        template<typename URBG>
        inline result_type operator()(URBG &rng) {
            return static_cast<result_type>(std::floor(m_d(rng)));
        }
#else
        template<typename URBG>
        inline result_type operator()(URBG &rng) {
            return m_id(rng);
        }
#endif

    private:
#ifdef __apple_build_version__
        const double m_a, m_b;
        std::uniform_real_distribution<double> m_d;
#else
        std::uniform_int_distribution<IntType> m_id;
#endif
    };
}// namespace distributions

#endif//ESERCIZI_LSN_UNIFORM_INT_HPP
