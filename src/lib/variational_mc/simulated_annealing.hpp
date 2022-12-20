//
// Created by Davide Nicoli on 24/10/22.
//

#ifndef ESERCIZI_LSN_SIMULATED_ANNEALING_HPP
#define ESERCIZI_LSN_SIMULATED_ANNEALING_HPP

#include <cstddef>
#include <iostream>
#include <tuple>

#include <indicators/progress_bar.hpp>

#include "samplers/MCMC/metropolis.hpp"

using namespace samplers::mcmc;
using namespace indicators;

namespace variational_mc {

    template<typename field>
    class [[maybe_unused]] LinearScheduler {
    public:
        LinearScheduler(field start, field end, size_t n_steps)
            : m_start(start), m_delta(end - start), m_n_steps(static_cast<field>(n_steps)),
              m_n_steps_int(m_n_steps) {}

        constexpr inline field value(size_t step) const {
            return m_start + m_delta * static_cast<field>(step) / m_n_steps;
        }

        [[nodiscard]] constexpr inline bool end(size_t step) const { return step == m_n_steps_int; }

        [[nodiscard]] ProgressBar pbar(std::string && /*var*/) const {
            return ProgressBar(option::ProgressType{ProgressType::decremental},
                               option::MaxProgress{m_n_steps_int});
        }

    private:
        const field m_start, m_delta, m_n_steps;
        const size_t m_n_steps_int;
    };


    template<typename field>
    class LogScheduler {
    public:
        LogScheduler(field start, field end, size_t n_steps)
            : m_start(start), m_frac(end / start), m_n_steps_1(static_cast<field>(n_steps - 1)),
              m_n_steps_int(n_steps) {
            if (n_steps < 1) throw std::runtime_error("Not enough steps");
        }

        constexpr inline field value(size_t step) const {
            return m_start * std::pow(m_frac, static_cast<field>(step) / m_n_steps_1);
        }

        [[nodiscard]] constexpr inline bool end(size_t step) const { return step >= m_n_steps_int; }


        [[nodiscard]] ProgressBar pbar(std::string &&var) const {
            return ProgressBar(option::ProgressType{ProgressType::decremental},
                               option::MaxProgress{m_n_steps_int},
                               option::PostfixText{"Log(" + var + ")"}, option::Lead("="));
        }

    private:
        const field m_start, m_frac, m_n_steps_1;
        const size_t m_n_steps_int;
    };

    // TODO: change the "field" param. It is reaaaally bad practice
    template<class StochasticLoss, class ConditionalTransition, typename field = double>
    class SimulatedAnnealing {
        using params_space = typename ConditionalTransition::StateSpace;

    public:
        SimulatedAnnealing(StochasticLoss &&loss_fn, ConditionalTransition q)
            : m_loss(std::forward<StochasticLoss>(loss_fn)), m_q(q) {}

        template<class DecayScheduler, typename ParamsIt, typename EnergyIt, typename TemperatureIt,
                 class URBG>
        void anneal(const params_space &p0, size_t explore_steps, const DecayScheduler &T_scheduler,
                    ParamsIt first_params, EnergyIt first_energy, EnergyIt first_uncert,
                    TemperatureIt first_T, URBG &rng) {
            ProgressBar pbar = T_scheduler.pbar("Temperature");
            pbar.print_progress();
            *first_params++ = p0;
            std::tie(*first_energy++, *first_uncert++) = m_loss(p0);
            size_t T_step = 0;
            const auto T0 = T_scheduler.value(0);
            *first_T++ = T0;
            typename EnergyIt::value_type benergy, berror;
            while (!T_scheduler.end(T_step)) {
                const auto T = T_scheduler.value(T_step);
                auto params_sampler = make_params_sampler(*std::prev(first_params), T);
                for (size_t step = 0; step < explore_steps; step++) {
                    std::tie(std::ignore, *first_params++, benergy, berror) =
                            params_sampler.step_p(rng);
                    *first_energy++ = -benergy * T;
                    *first_uncert++ = -berror * T;
                    *first_T++ = T;
                }
                T_step++;
                pbar.tick();
            }
        }

    private:
        StochasticLoss m_loss;
        ConditionalTransition m_q;

        class pdf {
        public:
            typedef typename ConditionalTransition::ProbSpace prob_space;

            explicit pdf(field T, StochasticLoss &loss) : m_T(T), m_loss(loss) {}

            constexpr inline auto logp(const params_space &x) {
                const auto [loss, uncert] = m_loss(x);
                return std::make_tuple(-static_cast<prob_space>(loss) / m_T,
                                       -static_cast<prob_space>(uncert) / m_T);
            }

        private:
            const prob_space m_T;
            StochasticLoss &m_loss;
        };

        inline auto make_params_sampler(const params_space &p0, field temperature) {
            return SAMetropolis(p0, pdf(temperature, m_loss), m_q);
        }
    };

}// namespace variational_mc

#endif//ESERCIZI_LSN_SIMULATED_ANNEALING_HPP
