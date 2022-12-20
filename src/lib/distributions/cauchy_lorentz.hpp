//
// Created by Davide Nicoli on 17/03/22.
//

#ifndef ESERCIZI_LSN_CAUCHY_LORENTZ_HPP
#define ESERCIZI_LSN_CAUCHY_LORENTZ_HPP

#include <cmath>
#include <random>

namespace distributions {
    template<typename value>
    class cauchy_lorentz {
    public:
        using result_type = value;
        using param_type = value;

        explicit cauchy_lorentz(const param_type mu, const param_type gamma)
            : m_mu(mu), m_gamma(gamma) {}

        template<class URBG>
        result_type operator()(URBG &rng) {
            return m_mu + m_gamma * tan((m_rd(rng) - .5) * M_PI);
        }

    private:
        param_type m_mu;
        param_type m_gamma;
        std::uniform_real_distribution<value> m_rd{0, 1};
    };

}// namespace distributions

#endif// ESERCIZI_LSN_CAUCHY_LORENTZ_HPP
