//
// Created by Davide Nicoli on 19/10/22.
//

#ifndef ESERCIZI_LSN_EXERCISES_HPP
#define ESERCIZI_LSN_EXERCISES_HPP

#include <cmath>
#include <valarray>

namespace distributions::ex08 {
    template<typename field>
    class Trial {
    public:
        typedef field StateSpace;
        typedef field prob_space;

        Trial(field mu, field sigma) : m_mu{mu}, m_sigma{sigma} {}

        /**
         * log |ψ(x)|^2
         * @param x
         * @return
         */
        constexpr prob_space logp(const field &x) {
            // This is a good way to deal with infinites...
            const auto xp = std::pow((x + m_mu) / m_sigma, 2);
            const auto xm = std::pow((x - m_mu) / m_sigma, 2);
            return 2 * std::log(std::exp(-xp / 2) + std::exp(-xm / 2));
            // ... this is not, as a1 can easily go to -inf and a2 to both +inf and -inf (sigma<<1)
            // const prob_space a1 = -(x * x + m_m2) / m_s2;
            // const prob_space a2 = std::log(1 + std::cosh(m_k1_ * x));
            // return std::log(2) + a1 + a2;
        }

    private:
        field m_mu, m_sigma;
    };
}// namespace distributions::ex08

namespace functions::ex08 {
    template<typename field>
    class Integrand {
    public:
        typedef field return_type;
        Integrand(field mu, field sigma) : m_mu(mu), m_sigma(sigma) {}

        /**
         * Hψ/ψ
         * @param x x
         */
        constexpr return_type operator()(const field &x) {
            const auto x2 = x * x, xmu = x * m_mu;
            return m_a1 - (x2 + m_m2 - 2 * xmu * std::tanh(xmu / m_s2)) / (2 * m_s4) + x2 * x2 -
                   2.5 * x2;
        }

    private:
        field m_mu, m_sigma, m_m2{m_mu * m_mu}, m_s2{m_sigma * m_sigma}, m_s4{m_s2 * m_s2},
                m_a1{1 / (2 * m_s2)};
    };
}// namespace functions::ex08

#endif//ESERCIZI_LSN_EXERCISES_HPP
