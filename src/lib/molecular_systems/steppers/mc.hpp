//
// Created by Davide Nicoli on 19/09/22.
//

#ifndef ESERCIZI_LSN_MS_STEPPER_MC_HPP
#define ESERCIZI_LSN_MS_STEPPER_MC_HPP

#include <cmath>
#include <random>
#include <vector>

#include "../algos.hpp"
#include "distributions/uniform_int.hpp"
#include "molecular_systems/system.hpp"

namespace molecular_systems::steppers {
    /**
     * Explores an LJMono system's states using the Metropolis-Hastings algorithm.
     * @tparam field System's numeric field.
     * @tparam tail_corrections Whether tail corrections to the thermodynamic variables are used.
     */
    template<typename field, bool tail_corrections, class URBG>
    class MC {
    public:
        typedef LJMono<field, tail_corrections, Ensamble::NVT> System;
        typedef field Field;
        MC(MC &) = delete;
        MC(const MC &) = delete;
        MC(MC &&) noexcept = default;

        /**
         * Constructor
         * @param n_particles Number of particles in the system.
         * @param displacement_diameter Diameter in which a tentative particle will be sampled.
         * @param rng Pointer to a random number generator.
         */
        explicit MC(size_t n_particles, field displacement_diameter, std::shared_ptr<URBG> rng)
            : m_n_particles(n_particles), m_particle(0, n_particles - 1),
              m_displacement(-displacement_diameter / 2, displacement_diameter / 2), m_rng(rng) {}

        /**
         * Samples a new state for the system using the Metropolis-Hastings algorithm.
         * @param system An LJMono system.
         */
        void step(System &system) {
            for (size_t i = 0; i < m_n_particles; i++) {
                // Sampling a particle to displace
                const size_t particle = m_particle(*m_rng);
                const auto old_position = system.m_positions[particle];
                // Sampled particle's energy
                const auto e_old = LJ_potential<tail_corrections>(
                        particle, old_position, system.m_positions, system.m_simulation);
                // Sampling a new position
                const auto newx = system.pbc(old_position[0] + m_displacement(*m_rng));
                const auto newy = system.pbc(old_position[1] + m_displacement(*m_rng));
                const auto newz = system.pbc(old_position[2] + m_displacement(*m_rng));
                // Displaced particle's energy
                const auto e_new = LJ_potential<tail_corrections>(
                        particle, {newx, newy, newz}, system.m_positions, system.m_simulation);
                // Metropolis' threshold
                const auto p = std::exp((e_old - e_new) / system.m_thermo.temperature);
                if (m_unit(*m_rng) <= p) {
                    //                    system.m_prev_positions.e_i[particle] = old_position[0];
                    //                    system.m_prev_positions.e_j[particle] = old_position[1];
                    //                    system.m_prev_positions.e_k[particle] = old_position[2];
                    system.m_positions.e_i[particle] = newx;
                    system.m_positions.e_j[particle] = newy;
                    system.m_positions.e_k[particle] = newz;
                    m_accepted_steps++;
                }
            }
            m_total_steps += m_n_particles;
            system.time_step();
        }

        // Acceptance rate during the last round
        [[nodiscard]] double acceptance_rate() const {
            return double(m_accepted_steps) / double(m_total_steps);
        }

        static constexpr bool tails = tail_corrections;
        static constexpr bool compute_forces = false;

    private:
        const size_t m_n_particles;
        distributions::uniform_int<size_t> m_particle;
        std::uniform_real_distribution<field> m_displacement;
        std::uniform_real_distribution<field> m_unit{0, 1};
        size_t m_accepted_steps{0};
        size_t m_total_steps{0};
        std::shared_ptr<URBG> m_rng;
    };
}// namespace molecular_systems::steppers

#endif//ESERCIZI_LSN_MS_STEPPER_MC_HPP
