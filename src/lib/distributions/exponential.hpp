//
// Created by Davide Nicoli on 17/03/22.
//

#ifndef ESERCIZI_LSN_EXPONENTIAL_HPP
#define ESERCIZI_LSN_EXPONENTIAL_HPP

#include <cmath>
#include <random>

namespace distributions {
    /**
     * An exponential distribution with an interface similar to C++'s <random> distributions
     */
    template<typename value>
    class exponential {
    public:
        typedef value result_type;
        typedef value param_type;

        explicit exponential(const param_type lambda) : m_lambda(lambda) {}

        /**
         * Produces exponentially distributed numbers.
         * @param rng A random number generator respecting C++'s UniformRandomBitGenerator requirements (https://en.cppreference.com/w/cpp/named_req/UniformRandomBitGenerator)
         */
        template<class URBG>
        result_type operator()(URBG &rng) {
            return -log(1 - m_rd(rng)) / m_lambda;
        }

    private:
        param_type m_lambda;
        std::uniform_real_distribution<value> m_rd{0, 1};
    };

}// namespace distributions

#endif// ESERCIZI_LSN_EXPONENTIAL_HPP
