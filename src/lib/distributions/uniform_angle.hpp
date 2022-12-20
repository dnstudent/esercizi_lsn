//
// Created by Davide Nicoli on 23/03/22.
//

#ifndef ESERCIZI_LSN_UNIFORM_ANGLE_HPP
#define ESERCIZI_LSN_UNIFORM_ANGLE_HPP

#include "goniom.hpp"

#include <cmath>
#include <random>

namespace distributions {
    /**
     * Samples uniformly transform direction so that the area density of generated points is uniform.
     * @tparam DirectionContainer The container in which to store sampled coordinates.
     */
    template<class DirectionContainer>
    class Uniform3DDirection {
    public:
        typedef DirectionContainer result_type;
        typedef typename result_type::value_type value_type;

        /**
         * Sampling method. Uses transform sin distribution on the polar angle to achieve surface density uniformity.
         * @param rng The random generator used to sample points.
         * @return A triplet of cartesian coordinates representing transform versor.
         */
        template<class URBG>
        result_type operator()(URBG &rng) {
            const auto theta = m_thd(rng);
            const auto sintheta = sin(theta);
            const auto phi = m_phid(rng);
            return {sintheta * cos(phi), sintheta * sin(phi), cos(theta)};
        }

    private:
        distributions::Sin<value_type> m_thd{};
        std::uniform_real_distribution<value_type> m_phid{0, 2 * M_PI};
    };
}// namespace distributions

#endif// ESERCIZI_LSN_UNIFORM_ANGLE_HPP
