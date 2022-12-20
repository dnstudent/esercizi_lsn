//
// Created by Davide Nicoli on 25/08/22.
//

#ifndef ESERCIZI_LSN_GAUSS_HPP
#define ESERCIZI_LSN_GAUSS_HPP

#include <cmath>
#include <random>
#include <valarray>

#include "transition.hpp"
#include "utils.hpp"

namespace transitions {
    template<typename prob_space, typename state_space>
    class GaussNear {
        using real_space = typename state_space::value_type;

    public:
        typedef state_space StateSpace;
        typedef prob_space ProbSpace;
        typedef void Symmetric;

        GaussNear(real_space stdev, size_t n_dimensions)
            : m_stdev{stdev}, m_ndimensions(n_dimensions) {}

        template<class URBG>
        inline StateSpace sample(const StateSpace &from, URBG &rng) {
            StateSpace to(m_ndimensions);
            std::transform(std::begin(from), std::end(from), std::begin(to),
                           [&](const auto xi) { return xi + m_offset(rng); });
            return to;
        }

        inline ProbSpace logp(const StateSpace &to, const StateSpace &from) const {
            const auto delta = to - from;
            return m_prefix - (delta * delta).sum() / m_2var;
        }

    protected:
        const real_space m_stdev, m_2var{2 * m_stdev * m_stdev};
        const size_t m_ndimensions;
        const real_space m_prefix{-real_space(m_ndimensions) *
                                  (0.5 * std::log(2 * M_PI) + std::log(m_stdev))};
        std::normal_distribution<real_space> m_offset{0, m_stdev};
    };

    template<typename prob_space, typename real_space, size_t N>
    class GaussNear<prob_space, std::array<real_space, N>> {
    public:
        typedef std::array<real_space, N> StateSpace;
        typedef prob_space ProbSpace;
        typedef void Symmetric;

        explicit GaussNear(real_space stdev) : m_stdev{stdev} {}

        template<class URBG>
        inline StateSpace sample(const StateSpace &from, URBG &rng) {
            StateSpace to{};
            std::transform(std::begin(from), std::end(from), std::begin(to),
                           [&](const auto xi) { return xi + m_offset(rng); });
            return to;
        }

        inline ProbSpace logp(const StateSpace &to, const StateSpace &from) const {
            const auto delta = to - from;
            return m_prefix - (delta * delta).sum() / m_2var;
        }

    protected:
        const real_space m_stdev, m_2var{2 * m_stdev * m_stdev};
        const real_space m_prefix{-real_space(N) * (0.5 * std::log(2 * M_PI) + std::log(m_stdev))};
        std::normal_distribution<real_space> m_offset{0, m_stdev};
    };

    template<typename prob_space>
    class GaussNear<prob_space, double> {
        using real_space = double;

    public:
        typedef double StateSpace;
        typedef prob_space ProbSpace;
        typedef void Symmetric;

        explicit GaussNear(real_space stdev) : m_stdev{stdev} {}

        template<class URBG>
        inline StateSpace sample(const StateSpace &from, URBG &rng) {
            return from + m_offset(rng);
        }

        inline ProbSpace logp(const StateSpace &to, const StateSpace &from) const {
            const auto delta = to - from;
            return m_prefix - (delta * delta) / m_2var;
        }

    protected:
        const real_space m_stdev, m_2var{2 * m_stdev * m_stdev};
        const real_space m_prefix{-(0.5 * std::log(2 * M_PI) + std::log(m_stdev))};
        std::normal_distribution<real_space> m_offset{0, m_stdev};
    };
}// namespace transitions

#endif//ESERCIZI_LSN_GAUSS_HPP
