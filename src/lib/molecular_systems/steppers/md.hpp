//
// Created by Davide Nicoli on 19/09/22.
//

#ifndef ESERCIZI_LSN_MS_STEPPER_MD_HPP
#define ESERCIZI_LSN_MS_STEPPER_MD_HPP

#include <memory>

#include "../algos.hpp"
#include "../system.hpp"
#include "collectors.hpp"

namespace molecular_systems::steppers {

    /**
     * Stepper used in exercise 04
     * @tparam field
     * @tparam tail_corrections
     */
    template<typename field, bool tail_corrections>
    class MD : public Stepper<LJMono<field, tail_corrections, Ensamble::NVE>,
                              MD<field, tail_corrections>> {
        typedef Vectors<field> V;
        using System = LJMono<field, tail_corrections, Ensamble::NVE>;
        using Base = Stepper<LJMono<field, tail_corrections, Ensamble::NVE>,
                             MD<field, tail_corrections>>;

    public:
        explicit MD(std::shared_ptr<System> system) : Base(system) {}

        void step() {
            Vectors<field> next_positions(this->m_system->m_simulation.n_particles);
            verlet_next_positions(this->m_system->m_positions, this->m_system->m_prev_positions,
                                  next_positions, this->m_system->m_forces,
                                  this->m_system->m_simulation.delta2,
                                  this->m_system->m_simulation.box_edge);
            verlet_next_velocities(
                    next_positions, this->m_system->m_prev_positions, this->m_system->m_velocities,
                    this->m_system->m_simulation.dbldelta, this->m_system->m_simulation.box_edge);
            this->m_system->m_prev_positions = std::move(this->m_system->m_positions);
            this->m_system->m_positions = std::move(next_positions);
            this->m_system->time_step();
        }
    };

    /**
     * Stepper used in exercise 07
     * @tparam field
     * @tparam tail_corrections
     */
    template<typename field, bool tail_corrections>
    class MD2 {
    public:
        typedef LJMono<field, tail_corrections, Ensamble::NVE> System;
        typedef field Field;
        MD2(MD2 &) = delete;
        MD2(const MD2 &) = delete;
        MD2(MD2 &&) noexcept = default;

        MD2() = default;

        /**
         * Evolves the system's state using the Verlet integration scheme
         * Needs to be followed by a call to "measures", as the software design
         * is very bad
         * @param system An LJMono system.
         */
        void step(System &system) {
            Vectors<field> next_positions(system.m_simulation.n_particles);
            verlet_next_positions(system.m_positions, system.m_prev_positions, next_positions,
                                  system.m_forces, system.m_simulation.delta2,
                                  system.m_simulation.box_edge);
            verlet_next_velocities(next_positions, system.m_prev_positions, system.m_velocities,
                                   system.m_simulation.dbldelta, system.m_simulation.box_edge);
            system.m_prev_positions = std::move(system.m_positions);
            system.m_positions = std::move(next_positions);
            system.time_step();
        }

        static constexpr bool tails = tail_corrections;
        static constexpr bool compute_forces = true;
    };
}// namespace molecular_systems::steppers

#endif//ESERCIZI_LSN_MS_STEPPER_MD_HPP
