//
// Created by Davide Nicoli on 03/05/22.
//

#ifndef ESERCIZI_LSN_GIBBS_HPP
#define ESERCIZI_LSN_GIBBS_HPP

#include <cmath>
#include <memory>
#include <random>

#include "models/ising/1D/ising.hpp"

using namespace ising;

namespace samplers::mcmc {

    /**
     * Gibbs MCMC sampler for the 1D Ising system. It shares a common interface with the SystemMetropolis sampler, so that the two of them are interchangeable.
     * @tparam var_space Ising model's thermo variables space.
     */
    template<typename var_space>
    class SystemGibbs {
        typedef typename Ising<var_space, 1>::StateSpace StateSpace;

    public:
        /**
         * Initializer.
         * @param system Pointer to an Ising system defined somewhere else.
         */
        explicit SystemGibbs(std::shared_ptr<Ising<var_space, 1>> system) : m_system(system) {}

        /**
         * Performs a Gibbs MC step. The result is stored in the system.
         * @param rng Random number generator.
         * @return Whether the proposed step was accepted or not (it is always accepted in this scheme).
         */
        template<class URBG>
        bool step(URBG &rng) {
            for (size_t i = 0; i < m_system->n_spins(); i++) {
                if (m_unif(rng) < p_k1(i)) m_system->set(i, spins::up);
                else
                    m_system->set(i, spins::down);
            }
            return m_system->n_spins();
        }


        template<class URBG>
        inline void warmup(size_t steps, URBG &rng) {
            for (size_t i = 0; i < steps; i++) { step(rng); }
        }

        /**
        * Performs a number of steps, each followed by an action.
        * @tparam Action An invocable class.
        * @param n_steps Number of steps.
        * @param action Action to perform; must take no arguments.
        * @param rng Random number generator.
        * @return Acceptance rate.
        */
        template<class Action, class URBG>
        double process(size_t n_steps, Action action, URBG &rng) {
            for (size_t i = 0; i < n_steps; i++) {
                step(rng);
                action();
            }
            return 1;
        }

    private:
        std::shared_ptr<Ising<var_space, 1>> m_system;
        std::uniform_real_distribution<double> m_unif{0, 1};

        [[nodiscard]] inline double p_k1(size_t k) const noexcept {
            const auto e_d = std::exp(-m_system->beta() * m_system->flip_dE(k));
            const auto num = m_system->state()[k] ? double(1) : e_d;
            return num / (double(1) + e_d);
        }
    };
}// namespace samplers::mcmc

#endif// ESERCIZI_LSN_GIBBS_HPP
