//
// Created by Davide Nicoli on 29/08/22.
//

#ifndef ESERCIZI_LSN_ISING_1D_SIMULATION_HPP
#define ESERCIZI_LSN_ISING_1D_SIMULATION_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <string_view>
#include <tuple>

#include <rapidcsv.h>

#include "estimators/estimators.hpp"
#include "estimators/mean.hpp"
#include "estimators/variance.hpp"
#include "ising.hpp"
#include "utils.hpp"
#include "variables.hpp"

/**
 * Classes used to perform runs. Originally I wanted to make them for generic models,
 * but I ended up specializing on the Ising 1D model. Some parts of the code may seem
 * a bit convoluted.
 */

using namespace estimators;

namespace fs = std::filesystem;
namespace csv = rapidcsv;

namespace ising::D1 {

    /**
     * Performs equilibration. It differs from the simulator in that it does not carry out block measurements
     * of thermodynamical variables, but only registers instantaneous values of proxy variables. It is
     * aimed at assessing autocorrelation times.
     * @tparam compute_H Compile time flag controlling whether the values of H will be computed and registered.
     * @tparam compute_s Compile time flag controlling whether the values of Sum_s will be computed and registered.
     * @tparam compute_s2 Compile time flag controlling whether the values of Sum_s2 will be computed and registered.
     * @tparam var_space The numerical field of the variables.
     * @tparam SystemSampler An MCMC sampler class specific for Ising 1D systems.
     */
    template<bool compute_H, bool compute_s, bool compute_s2, typename var_space,
             class SystemSampler>
    class Equilibrator {
    public:
        /**
         * Initializes the Equilibrator.
         * @param n_steps Number of steps for which proxy variables will be computed and stored.
         * @param model Ising model that will be equilibrated.
         */
        Equilibrator(size_t n_steps, std::shared_ptr<Ising<var_space, 1>> model)
            : m_n_steps(n_steps), m_model(model) {}

        /**
         * Performs a number of steps without computing proxy variables.
         * @param warmup_steps Number of steps.
         * @param rng Random number generator.
         */
        template<class URBG>
        inline void warmup(size_t warmup_steps, URBG &rng) {
            for (size_t i = 0; i < warmup_steps; i++) { m_sampler.step(rng); };
        }

        /**
         * Performs an equilibration run.
         * @param warmup_steps Number of steps at the beginning of the run for which no variable is computed.
         * @param rng Random number generator.
         */
        template<class URBG>
        void run(size_t warmup_steps, URBG &rng) {
            warmup(warmup_steps, rng);
            for (size_t sample = 0; sample < m_n_steps; sample++) {
                // Using compile-time switches for performance reasons.
                if constexpr (compute_H) m_vars.h[sample] = m_model->energy();
                if constexpr (compute_s) m_vars.sum_s[sample] = spin_sum(*m_model);
                if constexpr (compute_s2) m_vars.sum_s2[sample] = spin_sum2(*m_model);
                m_sampler.step(rng);
            }
        }

        inline void save_results(const fs::path &output_path) const {
            m_vars.save_data(output_path);
        }

        inline void save_state(const fs::path &output_path) const {
            m_model->save_state(output_path);
        }

    private:
        size_t m_n_steps;
        std::shared_ptr<Ising<var_space, 1>> m_model;
        SystemSampler m_sampler{m_model};
        block_proxy_vars<var_space> m_vars{m_n_steps};
    };

    /**
     * Performs simulations. It samples states of the Ising model through a given MCMC algorithm
     * and computes the estimations (with error) of specified thermodynamical variables performing
     * block statistics on instantaneous proxy data. Chosen proxies are the energy (h), the sum of spins
     * (Sum_s) and the sum of the product of spins (Sum_s2), as the th. variables of interest in
     * exercise 06 are defined on those. They could be abstracted away as it has been done with the
     * th. variables themselves, but this has a certain amount of complications.
     * @tparam H Whether the values of H need to be computed.
     * @tparam S Whether the values of Sum_s need to be computed.
     * @tparam S2 Whether the values of Sum_s2 need to be computed.
     * @tparam var_space The numerical field of the variables.
     * @tparam SystemSampler An MCMC sampler class specific for Ising 1D systems.
     * @tparam Ising1DThermoVars 1D Ising thermodynamical variables' classes.
     */
    template<bool H, bool S, bool S2, typename var_space, class SystemSampler,
             class... Ising1DThermoVars>
    class Simulator {
    public:
        /**
         * Initializer. Block size is fixed for performance reasons.
         * @param block_size Size of the blocks used to estimate the values of thermo variables.
         * @param model 1D ising model which will be evolved.
         * @param vars Sequence of variables to be estimated.
         */
        Simulator(size_t block_size, std::shared_ptr<Ising<var_space, 1>> model,
                  Ising1DThermoVars &&...vars)
            : m_block_size(block_size), m_model(model),
              m_thermovars(std::make_tuple(std::forward<Ising1DThermoVars>(vars)...)),
              m_cachevars(block_size) {}

        /**
         * Performs a series of MC steps and stores variables' estimates.
         * @param rng Random number generator.
         * @return Step acceptance rate: how many of the proposed steps are accepted?
         */
        template<class URBG>
        double block_estimates(URBG &rng) {
            const double acc_rate = m_sampler.process(
                    m_block_size,
                    [&]() {
                        measure(m_step % m_block_size);
                        m_step++;
                    },
                    rng);
            //for (size_t i = 0; i < m_block_size; i++) { step(rng); }
            compute_store_estimates(m_cachevars, m_thermovars, m_thermo_outputs);
            return acc_rate;
        }

        /**
         * Performs a simulation.
         * @param n_blocks How many blocks should be processed.
         * @param n_warmup Number of steps to be used as warmup time.
         * @param rng Random number generator.
         */
        template<class URBG>
        void run(size_t n_blocks, size_t n_warmup, URBG &rng) {
            // for (size_t i = 0; i < n_warmup; i++) m_sampler.step(rng);
            m_sampler.warmup(n_warmup, rng);
            for (size_t block = 0; block < n_blocks; block++) block_estimates(rng);
        }

        void save_results(const fs::path &output_path) {
            csv::Document table;
            rec_store_results(table);
            if (!fs::exists(output_path.parent_path()))
                fs::create_directories(output_path.parent_path());
            table.RemoveColumn(table.GetColumnCount() - 1);
            table.Save(output_path);
        }

        void save_state(const fs::path &output_path) const { m_model->save_state(output_path); }


    private:
        size_t m_block_size, m_step{0};
        std::shared_ptr<Ising<var_space, 1>> m_model;
        SystemSampler m_sampler{m_model};
        std::tuple<Ising1DThermoVars...> m_thermovars;
        // Tuple of vectors of pairs
        std::array<std::pair<std::vector<var_space>, std::vector<var_space>>,
                   sizeof...(Ising1DThermoVars)>
                m_thermo_outputs{};
        block_proxy_vars<var_space> m_cachevars;


        /**
         * Measure proxy variables. Using compile-time flags for performance reasons.
         * @param sample Sample number.
         */
        void measure(size_t sample) {
            if constexpr (H) { m_cachevars.h[sample] = m_model->energy(); }
            if constexpr (S) {
                m_cachevars.sum_s[sample] = spin_sum(*m_model);
                if constexpr (S2)
                    m_cachevars.sum_s2[sample] =
                            m_cachevars.sum_s[sample] * m_cachevars.sum_s[sample];
            }
            if constexpr (S2 && !S) {
                // ;
                m_cachevars.sum_s2[sample] = spin_sum2(*m_model);
            }
        }

        /**
         * Utility to store estimations and errors in a rapidcsv table.
         * @tparam I Internal use.
         * @param table Rapidcsv table.
         */
        template<size_t I = 0>
        inline constexpr void rec_store_results(csv::Document &table) {
            if constexpr (I == std::tuple_size_v<decltype(m_thermovars)>) return;
            else {
                const auto var_name = std::get<I>(m_thermovars).name();
                table.InsertColumn(2 * I, std::get<I>(m_thermo_outputs).first,
                                   var_name + "_estimate");
                table.InsertColumn(2 * I + 1, std::get<I>(m_thermo_outputs).second,
                                   var_name + "_error");
                rec_store_results<I + 1>(table);
            }
        }
    };
}// namespace ising::D1

#endif//ESERCIZI_LSN_ISING_1D_SIMULATION_HPP
