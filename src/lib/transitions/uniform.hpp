//
// Created by Davide Nicoli on 23/08/22.
//

#ifndef ESERCIZI_LSN_UNIFORM_HPP
#define ESERCIZI_LSN_UNIFORM_HPP

#include <limits>
#include <random>
#include <valarray>

#include "transition.hpp"
#include "utils.hpp"

namespace transitions {
    //    template<typename prob_space, typename state_space>
    //    class UniformNear {};

    /// Transition using a box-volume uniform distribution
    /// \tparam prob_space Numeric field used for probabilities
    /// \tparam real_space Numeric field used as base for the vector space
    template<typename prob_space, typename real_space>
    class UniformNear {
    public:
        using StateSpace = real_space;
        using ProbSpace = prob_space;
        /// Symmetric flag; speeds up computation in the Metropolis-Hastings algorithm
        typedef void Symmetric;

        /// Constructor of the transition
        /// \param radius The uniform box half-edge
        /// \param n_dimensions The dimensionality of the vector space
        explicit UniformNear(real_space radius) : m_radius{radius} {}


        template<class URBG>
        inline StateSpace sample(const StateSpace &from, URBG &rng) {
            return from + m_offset(rng);
        }

        inline ProbSpace logp(const StateSpace &to, const StateSpace &from) const {
            if (std::abs(to - from) > m_radius) return -std::numeric_limits<ProbSpace>::infinity();
            return m_loginvnorm;
        }

    protected:
        const real_space m_radius;
        const ProbSpace m_invnorm{ProbSpace{1} / (2 * ProbSpace(m_radius))},
                m_loginvnorm{-std::log(2 * ProbSpace(m_radius))};
        std::uniform_real_distribution<real_space> m_offset{-m_radius, m_radius};
    };


    /// Transition using a box-volume uniform distribution
    /// \tparam prob_space Numeric field used for probabilities
    /// \tparam real_space Numeric field used as base for the vector space
    template<typename prob_space, typename real_space>
    class UniformNear<prob_space, std::valarray<real_space>> {
    public:
        typedef std::valarray<real_space> StateSpace;
        typedef prob_space ProbSpace;
        /// Symmetric flag; speeds up computation in the Metropolis-Hastings algorithm
        typedef void Symmetric;

        /// Constructor of the transition
        /// \param radius The uniform box half-edge
        /// \param n_dimensions The dimensionality of the vector space
        UniformNear(real_space radius, size_t n_dimensions)
            : m_radius{radius}, m_ndimensions(n_dimensions) {}


        template<class URBG>
        inline StateSpace sample(const StateSpace &from, URBG &rng) {
            StateSpace next(m_ndimensions);
            std::transform(std::begin(from), std::end(from), std::begin(next),
                           [&](const auto xi) { return xi + m_offset(rng); });
            return next;
        }

        inline ProbSpace logp(const StateSpace &to, const StateSpace &from) const {
            for (size_t i = 0; i < m_ndimensions; i++) {
                if (std::abs(to[i] - from[i]) > m_radius)
                    return -std::numeric_limits<ProbSpace>::infinity();
            }
            return m_loginvnorm;
        }

    protected:
        const real_space m_radius, m_radius_2{m_radius * m_radius};
        const size_t m_ndimensions;
        const ProbSpace m_invnorm{ProbSpace{1} / std::pow(2 * ProbSpace(m_radius), m_ndimensions)},
                m_loginvnorm{-ProbSpace(m_ndimensions) * std::log(2 * ProbSpace(m_radius))};
        std::uniform_real_distribution<real_space> m_offset{-m_radius, m_radius};
    };

    /// Transition using a box-volume uniform distribution
    /// \tparam prob_space Numeric field used for probabilities
    /// \tparam real_space Numeric field used as base for the vector space
    template<typename prob_space, typename real_space, size_t N>
    class UniformNear<prob_space, std::array<real_space, N>> {
    public:
        typedef std::valarray<real_space> StateSpace;
        typedef prob_space ProbSpace;
        /// Symmetric flag; speeds up computation in the Metropolis-Hastings algorithm
        typedef void Symmetric;

        /// Constructor of the transition
        /// \param radius The uniform box half-edge
        /// \param n_dimensions The dimensionality of the vector space
        explicit UniformNear(real_space radius) : m_radius{radius} {}


        template<class URBG>
        inline StateSpace sample(const StateSpace &from, URBG &rng) {
            StateSpace next;
            std::transform(std::begin(from), std::end(from), std::begin(next),
                           [&](const auto xi) { return xi + m_offset(rng); });
            return next;
        }

        inline ProbSpace logp(const StateSpace &to, const StateSpace &from) const {
            for (size_t i = 0; i < N; i++) {
                if (std::abs(to[i] - from[i]) > m_radius)
                    return -std::numeric_limits<ProbSpace>::infinity();
            }
            return m_loginvnorm;
        }

    protected:
        const real_space m_radius, m_radius_2{m_radius * m_radius};
        const ProbSpace m_invnorm{ProbSpace{1} / std::pow(2 * ProbSpace(m_radius), N)},
                m_loginvnorm{-ProbSpace(N) * std::log(2 * ProbSpace(m_radius))};
        std::uniform_real_distribution<real_space> m_offset{-m_radius, m_radius};
    };


    template<typename prob_space, typename real_space>
    class Uniform {
    public:
        typedef real_space StateSpace;
        typedef prob_space ProbSpace;
        /// Symmetric flag; speeds up computation in the Metropolis-Hastings algorithm
        typedef void Symmetric;

        /// Constructor of the transition
        /// \param radius The uniform box half-edge
        /// \param n_dimensions The dimensionality of the vector space
        Uniform(real_space a, real_space b) : m_a{a}, m_b(b) {}


        template<class URBG>
        inline StateSpace sample(const StateSpace & /*from*/, URBG &rng) {
            return m_sampler(rng);
        }

        inline ProbSpace logp(const StateSpace &to, const StateSpace & /*from*/) const {
            if (to < m_a || to > m_b) return -std::numeric_limits<ProbSpace>::infinity();
            return m_loginvnorm;
        }

    protected:
        const real_space m_a, m_b;
        const ProbSpace m_invnorm{ProbSpace{1} / static_cast<ProbSpace>(m_b - m_a)},
                m_loginvnorm{-static_cast<ProbSpace>(std::log(m_b - m_a))};
        std::uniform_real_distribution<real_space> m_sampler{m_a, m_b};
    };

}// namespace transitions


#endif//ESERCIZI_LSN_UNIFORM_HPP
