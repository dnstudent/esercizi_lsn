//
// Created by Davide Nicoli on 21/09/22.
//

#ifndef ESERCIZI_LSN_MS_STEPPER_HPP
#define ESERCIZI_LSN_MS_STEPPER_HPP

#include <array>
#include <filesystem>
#include <memory>
#include <optional>
#include <tuple>
#include <utility>

#include "../system.hpp"

namespace fs = std::filesystem;

namespace molecular_systems::steppers {

    template<class System, class Impl>
    class Stepper {
        using field = typename System::Field;
        typedef Vectors<field> V;

    public:
        explicit Stepper(std::shared_ptr<System> system) : m_system(system) {}


        /**
        * Computes mean and variance block estimates for each thermodynamical variable and checkpoints molecular positions in xyz format
        * @param save_every_n_frames Checkpoint interval
        * @param dir Directory where frames will be saved
        */
        auto block_estimates(size_t save_every, const std::optional<fs::path> &dir) {
            assert(save_every == 0 || dir.has_value());
            std::array<std::vector<field>, System::N_VARS()> block_data;
            for (auto &vec: block_data) { vec.resize(m_system->m_simulation.block_size); }
            auto stepper = static_cast<Impl &>(*this);
            for (size_t i = 0UL; i < m_system->m_simulation.block_size; i++) {
                m_system->measures(block_data[0][i], block_data[1][i], block_data[2][i],
                                   block_data[3][i], block_data[4][i]);
                stepper.step();
                if ((save_every > 0) && ((m_frame_counter % save_every) == 0)) {
                    m_system->save_xyz_positions(dir.value(), m_frame_counter);
                }
                m_frame_counter++;
            }
            std::array<std::tuple<field, field, field>, System::N_VARS()> results;
            for (size_t var = 0UL; var < System::N_VARS(); var++) {
                results[var] = m_estimators[var](block_data[var].cbegin(), block_data[var].cend());
            }
            return results;
        }


        /**
        * Computes mean and variance block estimates for each thermodynamical variable and checkpoints molecular positions in xyz format
        * @param save_every_n_frames Checkpoint interval
        * @param dir Directory where frames will be saved
        */
        template<Variable... vars>
        auto block_estimates2(size_t save_every, const std::optional<fs::path> &dir) {
            assert(save_every == 0 || dir.has_value());
            MeasureOutputs<field, vars...> outputs(m_system->m_simulation.block_size);
            std::array<std::vector<field>, System::N_VARS()> block_data;
            for (auto &vec: block_data) { vec.resize(m_system->m_simulation.block_size); }
            auto stepper = static_cast<Impl &>(*this);
            for (size_t i = 0UL; i < m_system->m_simulation.block_size; i++) {
                m_system->measures(block_data[0][i], block_data[1][i], block_data[2][i],
                                   block_data[3][i], block_data[4][i]);
                stepper.step();
                if ((save_every > 0) && ((m_frame_counter % save_every) == 0)) {
                    m_system->save_xyz_positions(dir.value(), m_frame_counter);
                }
                m_frame_counter++;
            }
            std::array<std::tuple<field, field, field>, System::N_VARS()> results;
            for (size_t var = 0UL; var < System::N_VARS(); var++) {
                results[var] = m_estimators[var](block_data[var].cbegin(), block_data[var].cend());
            }
            return results;
        }


        auto values_t() {
            auto &stepper = static_cast<Impl &>(*this);
            std::array<std::vector<field>, System::N_VARS()> measures;
            for (auto &vec: measures) { vec.resize(m_system->m_simulation.block_size); }
            for (size_t t = 0; t < m_system->m_simulation.block_size; t++) {
                m_system->measures(measures[0][t], measures[1][t], measures[2][t], measures[3][t],
                                   measures[4][t]);
                stepper.step();
            }
            return measures;
        }

        inline void checkpoint(const fs::path &positions_path,
                               const fs::path &velocities_path) const {
            m_system->save_configurations(positions_path, velocities_path);
        }

    protected:
        std::shared_ptr<System> m_system;
        size_t m_frame_counter{0};
        // Mean estimators for the thermodynamical variables
        std::array<SampleProgAvg<field>, System::N_VARS()> m_estimators{};
    };

    namespace detail {
        template<class C>
        struct output {
            typedef typename C::Output type;
        };

        template<class C>
        using output_t = typename output<C>::type;

        template<template<class...> class Extractor, class TupleLike>
        struct extract {};

        template<template<class...> class Extractor, template<class...> class TupleLike,
                 class... Elems>
        struct extract<Extractor, TupleLike<Elems...>> {
            typedef TupleLike<Extractor<Elems>...> type;
        };

        template<template<class...> class Extractor, class TupleLike>
        using extract_t = typename extract<Extractor, TupleLike>::type;


    }// namespace detail


    template<class Stepper, class Estimators, Variable... vars>
    class BlockStats {
        using System = typename Stepper::System;
        using field = typename Stepper::Field;
        using Outs = MeasureOutputs<field, vars...>;
        using EstimatorsOuts = detail::extract_t<detail::output_t, Estimators>;
        static_assert(std::tuple_size<Estimators>() == sizeof...(vars));

    public:
        BlockStats(Stepper &&stepper, Estimators &&estimators, size_t block_size)
            : m_stepper(std::forward<Stepper>(stepper)),
              m_estimators(std::forward<Estimators>(estimators)), m_measures(block_size) {}

        auto statistics(System &system) {
            for (size_t i = 0UL; i < system.m_simulation.block_size; i++) {
                system.template measures<Stepper::compute_forces>(m_measures);
                m_stepper.step(system);
            }
            EstimatorsOuts results;
            utils::tuple_transform(
                    [](const auto &var_measures, auto &var_estimator) {
                        return var_estimator(var_measures.cbegin(), var_measures.cend());
                    },
                    /*outs=*/results, /*ins=*/m_measures.all_measures(), m_estimators);
            m_measures.clear();
            return results;
        }

    private:
        Stepper m_stepper;
        Estimators m_estimators;
        Outs m_measures;
    };

    template<class Stepper, Variable... vars>
    class StepSampler {
        using System = typename Stepper::System;
        using field = typename Stepper::Field;
        using Outs = MeasureOutputs<field, vars...>;

    public:
        explicit StepSampler(Stepper &&stepper) : m_stepper(std::forward<Stepper>(stepper)) {}

        auto sample(System &system, size_t n_steps) {
            Outs measures(n_steps);
            for (size_t i = 0UL; i < n_steps; i++) {
                system.template measures<Stepper::compute_forces>(measures);
                m_stepper.step(system);
            }
            return measures;
        }

        Stepper m_stepper;
    };

}// namespace molecular_systems::steppers

#endif