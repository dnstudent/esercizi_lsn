//
// Created by Davide Nicoli on 17/03/22.
//

#ifndef ESERCIZI_LSN_GONIOM_HPP
#define ESERCIZI_LSN_GONIOM_HPP

#include <cmath>
#include <random>

namespace distributions {

    // p(x) = pi/2 Cos(x*pi/2) x in [0, 1]
    template<class RealType>
    class Cos {
    public:
        typedef RealType result_type;
        template<class URBG>
        constexpr result_type operator()(URBG &rng) const {
            return (2 / M_PI) * asin(m_dist(rng));
        }

    private:
        std::uniform_real_distribution<RealType> m_dist{0, 1};
    };

    // p(x) = sin(x) / 2; x in [0, pi]
    template<class RealType = double>
    class Sin {
    public:
        typedef RealType result_type;

        template<class URBG>
        constexpr result_type operator()(URBG &rng) {
            return acos(1 - 2 * m_dist(rng));
        }

    private:
        std::uniform_real_distribution<RealType> m_dist{0, 1};
    };
}// namespace distributions

#endif// ESERCIZI_LSN_GONIOM_HPP
