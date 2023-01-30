//
// Created by Davide Nicoli on 23/08/22.
//

#ifndef ESERCIZI_LSN_MCMC_METROPOLIS_HPP
#define ESERCIZI_LSN_MCMC_METROPOLIS_HPP

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <memory>
#include <random>
#include <tuple>
#include <type_traits>
#include <valarray>

#include "models/ising/1D/ising.hpp"
#include "transitions/transition.hpp"

using namespace ising;

/**
 * Structs to check at compile time whether a given operation has been defined as "Symmetric"
 * @tparam T Operation as a class/struct
 */
template<typename T, class = void>
struct is_symmetric : std::false_type {};

template<typename T>
struct is_symmetric<T, std::void_t<typename T::Symmetric>> : std::true_type {};


template<typename T, class = void>
struct is_stochastic : std::false_type {};

template<typename T>
struct is_stochastic<T, std::void_t<typename T::Stochastic>> : std::true_type {};

namespace samplers::mcmc {
    /**
     * Metropolis-Hastings sampler. It is quite generic.
     * @tparam state_space Space in which the algorithm will sample points.
     * @tparam PDF Probability density function: a struct/class which defines a "logp" (logarithm of a point's p. density value) method on points in state_space.
     * @tparam ConditionalTransition A struct/class that defines a "sample" (produce a new point from a given one) method on a point in state_space and a "logp" method on two points in state_space.
     */
    template<typename PDF, class ConditionalTransition>
    class Metropolis {
    public:
        typedef typename ConditionalTransition::StateSpace StateSpace;
        typedef typename PDF::prob_space ProbSpace;
        /**
         * Initializer.
         * @param start Point from which sampling begins.
         * @param pdf Probability density.
         * @param transition Conditional transition.
         */
        Metropolis(const StateSpace &start, PDF pdf, ConditionalTransition transition)
            : m_state(start), m_pdf(pdf), m_q(transition) {}

        /**
         * Performs a Metropolis MC step.
         * @param rng Random number generator.
         * @return Pair ("proposed step was accepted", "next point").
         */
        template<class URBG>
        std::tuple<bool, StateSpace> step(URBG &rng) {
            // Sampling a candidate point
            const StateSpace candidate = m_q.sample(m_state, rng);
            const auto candidate_logp = m_pdf.logp(candidate);
            if constexpr (is_stochastic<PDF>::value) { m_state_logp = m_pdf.logp(m_state); }
            ProbSpace step_log_prob;
            // Checking at compile time whether the transition has been marked as symmetric;
            // if it is unneded calculations are not performed.
            if constexpr (is_symmetric<ConditionalTransition>::value) {
                step_log_prob = candidate_logp - m_state_logp;
            } else {
                step_log_prob = candidate_logp + m_q.logp(m_state, candidate) - m_state_logp -
                                m_q.logp(candidate, m_state);
            }
            const bool accepted = m_uniform(rng) < std::exp(step_log_prob);
            if (accepted) {
                m_state = std::move(candidate);
                m_state_logp = std::move(candidate_logp);
            }
            return {accepted, m_state};
        }

        template<class URBG>
        std::tuple<bool, StateSpace, ProbSpace> step_p(URBG &rng) {
            // Sampling a candidate point
            const StateSpace candidate = m_q.sample(m_state, rng);
            const auto candidate_logp = m_pdf.logp(candidate);
            ProbSpace step_log_prob;
            // Checking at compile time whether the transition is has been marked as symmetric;
            // if it is unneded calculations are not performed.
            if constexpr (is_symmetric<ConditionalTransition>::value) {
                step_log_prob = candidate_logp - m_state_logp;
            } else {
                step_log_prob = candidate_logp + m_q.logp(m_state, candidate) - m_state_logp -
                                m_q.logp(candidate, m_state);
            }
            const bool accepted = m_uniform(rng) < std::exp(step_log_prob);
            if (accepted) {
                m_state = std::move(candidate);
                m_state_logp = std::move(candidate_logp);
            }
            return {accepted, m_state, m_state_logp};
        }

        template<class URBG>
        void warmup(size_t steps, URBG &rng) {
            for (size_t i = 0; i < steps; i++) { step(rng); }
        }

        /**
         * Performs a number of MC steps and stores the result.
         * @param first Beginning of output.
         * @param last Output's past-the-end iterator.
         * @param rng Random number generator.
         * @return Acceptance rate.
         */
        template<typename It, class URBG>
        double sample(It first, It last, URBG &rng) {
            std::generate(first, last, [&]() {
                auto [was_accepted, point] = step(rng);
                if (was_accepted) { m_accepted++; }
                return point;
            });
            m_processed += size_t(std::distance(first, last));
            return double(m_accepted) / double(m_processed);
        }

        template<typename PIt, typename EIt, class URBG>
        double sample_p(PIt first_point, PIt last_point, EIt first_prob, URBG &rng) {
            bool was_accepted;
            while (first_point != last_point) {
                std::tie(was_accepted, *first_point++, *first_prob++) = step_p(rng);
                if (was_accepted) m_accepted++;
            }
            m_processed += size_t(std::distance(first_point, last_point));
            return double(m_accepted) / double(m_processed);
        }

        /**
         * Performs a number of MC steps and stores the transformed result.
         * @tparam Transform Transformation from StateSpace to the value pointed to by It.
         * @param first Beginning of output.
         * @param last Output's past-the-end iterator.
         * @param rng Random number generator.
         * @param f Tranformation.
         * @return Acceptance rate.
         */
        template<typename It, class URBG, class Transform>
        double sample(It first, It last, URBG &rng, Transform f) {
            std::generate(first, last, [&]() {
                auto [was_accepted, point] = step(rng);
                if (was_accepted) { m_accepted++; }
                return f(point);
            });
            m_processed += size_t(std::distance(first, last));
            return double(m_accepted) / double(m_processed);
        }


    protected:
        StateSpace m_state;
        PDF m_pdf;
        ConditionalTransition m_q;
        // Persisting the state's logprob in case its computation is long
        ProbSpace m_state_logp{m_pdf.logp(m_state)};
        std::uniform_real_distribution<ProbSpace> m_uniform{0, 1};
        size_t m_accepted{0};
        size_t m_processed{0};
    };

    template<typename StochasticLoss, class ConditionalTransition>
    class SAMetropolis {
    public:
        typedef typename ConditionalTransition::StateSpace StateSpace;
        typedef typename StochasticLoss::prob_space ProbSpace;
        /**
         * Initializer.
         * @param start Point from which sampling begins.
         * @param loss Probability density.
         * @param transition Conditional transition.
         */
        SAMetropolis(const StateSpace &start, StochasticLoss loss, ConditionalTransition transition)
            : m_state(start), m_loss(loss), m_q(transition) {}


        template<class URBG>
        std::tuple<bool, StateSpace, ProbSpace, ProbSpace> step_p(URBG &rng) {
            // Sampling a candidate point
            const StateSpace candidate = m_q.sample(m_state, rng);
            const auto [candidate_logp, c_uncert] = m_loss.logp(candidate);
            auto [state_logp, s_uncert] = m_loss.logp(m_state);
            ProbSpace step_log_prob;
            // Checking at compile time whether the transition is has been marked as symmetric;
            // if it is unneded calculations are not performed.
            if constexpr (is_symmetric<ConditionalTransition>::value) {
                step_log_prob = candidate_logp - state_logp;
            } else {
                step_log_prob = candidate_logp + m_q.logp(m_state, candidate) - state_logp -
                                m_q.logp(candidate, m_state);
            }
            const bool accepted = m_uniform(rng) < std::exp(step_log_prob);
            if (accepted) {
                m_state = std::move(candidate);
                state_logp = candidate_logp;
                s_uncert = c_uncert;
            }
            return {accepted, m_state, state_logp, s_uncert};
        }

        template<class URBG>
        void warmup(size_t steps, URBG &rng) {
            for (size_t i = 0; i < steps; i++) { step_p(rng); }
        }


    protected:
        StateSpace m_state;
        StochasticLoss m_loss;
        ConditionalTransition m_q;
        // Persisting the state logp in case its computation is long
        std::uniform_real_distribution<ProbSpace> m_uniform{0, 1};
        size_t m_accepted{0};
        size_t m_processed{0};
    };

    /**
     * Metropolis-Hastings sampler on the 1D Ising model.
     * @tparam var_space Ising model's thermo variables space
     */
    template<typename var_space>
    class SystemMetropolis {
        using System = Ising<var_space, 1>;
        typedef typename System::StateSpace StateSpace;

    public:
        /**
         * Initializer.
         * @param system Pointer to a system declared somewhere else.
         */
        explicit SystemMetropolis(std::shared_ptr<System> system)
            : m_system(system), m_candidates(system->n_spins()) {
            std::iota(m_candidates.begin(), m_candidates.end(), size_t(0));
        }


        /**
         * Performs a Metropolis MC step. The result is stored in the system.
         * @param rng Random number generator.
         * @return The number of accepted steps.
         */
        template<class URBG>
        size_t step(URBG &rng) {
            std::shuffle(m_candidates.begin(), m_candidates.end(), rng);
            size_t n_accepted = 0;
            for (auto candidate: m_candidates) {
                const auto flip_log_prob = m_system->template flip_logp<double>(candidate);
                const bool accepted = m_uniform(rng) < std::exp(flip_log_prob);
                if (accepted) {
                    n_accepted++;
                    m_system->flip(candidate);
                }
            }
            return n_accepted;
        }

        template<class URBG>
        inline void warmup(size_t steps, URBG &rng) {
            for (size_t i = 0; i < steps; i++) step(rng);
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
                const auto n_accepted = step(rng);
                action();
                m_accepted += n_accepted;
            }
            m_processed += m_system->n_spins() * n_steps;
            return double(m_accepted) / double(m_processed);
        }

    private:
        std::shared_ptr<System> m_system;
        std::uniform_real_distribution<double> m_uniform{0, 1};
        size_t m_accepted{0};
        size_t m_processed{0};
        std::vector<size_t> m_candidates;
    };
}// namespace samplers::mcmc


#endif//ESERCIZI_LSN_MCMC_METROPOLIS_HPP
