//
// Created by Davide Nicoli on 24/07/22.
//

#ifndef ESERCIZI_LSN_MS_UTILS_HPP
#define ESERCIZI_LSN_MS_UTILS_HPP

#include <cmath>
#include <cstddef>
#include <random>
#include <valarray>

#include "data_types/settings.hpp"
#include "data_types/vectors.hpp"

/**
 * Periodic boundary conditions algorithm
 * @param x A coordinate component
 * @param box_edge The volume edge
 */
template<typename field>
inline field PBC(field x, field box_edge) {
    return x - box_edge * std::rint(x / box_edge);
}

/**
 * Computes previous positions from current positions and velocities, and stores it in transform given variable
 * @param positions Current positions
 * @param velocities Current velocities
 * @param dt Time step
 * @param box_edge Periodic boundary size
 * @param previous_positions_out Variables where the result will be stored
 */
template<typename field>
inline void compute_previous_positions(Vectors<field> &positions, Vectors<field> &velocities,
                                       field dt, field box_edge,
                                       Vectors<field> &previous_positions_out) {
    // convoluted approach to compute previous positions, but it avoids multiple allocations
    // without having to extend Vectors
    previous_positions_out = velocities;
    previous_positions_out *= dt;
    previous_positions_out -= positions;
    previous_positions_out.apply([&](auto &x) { x = PBC(-x, box_edge); });
}

/**
 * Generates particle velocities following transform Maxwell-Boltzmann distribution, using reduced units (k_b = 1, m = 1)
 * @param velocities Variable where the result will be stored
 * @param temperature The system temperature in reduced units
 * @param rng The random number generator to use
 */
template<typename field, class URBG>
void generate_velocities(Vectors<field> &velocities, field temperature, URBG &rng) {
    const auto n_particles = velocities.e_i.size();
    std::normal_distribution<field> gauss(0, std::sqrt(temperature));
    velocities.apply([&](auto &v) { v = gauss(rng); });
    velocities -= velocities.mean();
    const field v_scale_factor =
            std::sqrt(temperature / (velocities.full_norm2() / field(3 * n_particles)));
    velocities *= v_scale_factor;
}

/**
 * Computes next positions (on a single axis) using the Verlet integration scheme for transform list of particles
 * @param current_positions A single component of the particles positions
 * @param previous_positions A single component of the particles previous positions
 * @param next_positions The variable where results will be stored
 * @param forces A single component of the forces acting on each particle
 * @param dt2 The time step squared
 * @param box_edge The periodic boundary condition size
 */
template<typename field>
inline void verlet_next_positions(const std::valarray<field> &current_positions,
                                  const std::valarray<field> &previous_positions,
                                  std::valarray<field> &next_positions,
                                  const std::valarray<field> &forces, field dt2, field box_edge) {
    next_positions = 2 * current_positions - previous_positions + forces * dt2;
    std::for_each(std::begin(next_positions), std::end(next_positions),
                  [=](auto &x) { x = PBC(x, box_edge); });
}

/**
 * Computes next positions (on three axes) using the Verlet integration scheme for transform list of particles
 * @param current_positions The particles positions
 * @param previous_positions The particles previous positions
 * @param next_positions The variable where results will be stored
 * @param forces The forces acting on each particle
 * @param dt2 The time step squared
 * @param box_edge The periodic boundary condition size
 */
template<typename field>
inline void verlet_next_positions(Vectors<field> &current_positions,
                                  Vectors<field> &previous_positions,
                                  Vectors<field> &next_positions, const Vectors<field> &forces,
                                  field dt2, field box_edge) {
    verlet_next_positions(current_positions.e_i, previous_positions.e_i, next_positions.e_i,
                          forces.e_i, dt2, box_edge);
    verlet_next_positions(current_positions.e_j, previous_positions.e_j, next_positions.e_j,
                          forces.e_j, dt2, box_edge);
    verlet_next_positions(current_positions.e_k, previous_positions.e_k, next_positions.e_k,
                          forces.e_k, dt2, box_edge);
}

/**
* Computes the velocities (on a single axis) of transform list of particles using the Verlet integration scheme
* @param next_positions The particles next positions
* @param previous_positions The particles previous positions
* @param next_velocities The variable where results will be stored
* @param dbldt Double the time step
* @param box_edge The periodic boundary condition size
 */
template<typename field>
inline void verlet_next_velocities(const std::valarray<field> &next_positions,
                                   const std::valarray<field> &previous_positions,
                                   std::valarray<field> &next_velocities, field dbldt,
                                   field box_edge) {
    next_velocities = (next_positions - previous_positions);
    std::for_each(std::begin(next_velocities), std::end(next_velocities),
                  [=](field &d) { d = PBC(d, box_edge) / dbldt; });
}

/**
* Computes the velocities (on three axes) of transform list of particles using the Verlet integration scheme
* @param next_positions The particles next positions
* @param previous_positions The particles previous positions
* @param next_velocities The variable where results will be stored
* @param dbldt Double the time step
* @param box_edge The periodic boundary condition size
 */
template<typename field>
inline void verlet_next_velocities(const Vectors<field> &next_positions,
                                   const Vectors<field> &previous_positions,
                                   Vectors<field> &next_velocities, field dbldt, field box_edge) {
    verlet_next_velocities(next_positions.e_i, previous_positions.e_i, next_velocities.e_i, dbldt,
                           box_edge);
    verlet_next_velocities(next_positions.e_j, previous_positions.e_j, next_velocities.e_j, dbldt,
                           box_edge);
    verlet_next_velocities(next_positions.e_k, previous_positions.e_k, next_velocities.e_k, dbldt,
                           box_edge);
}


/**
 * Computes the Lennard-Jones potential (in reduced units) acted by particles in "positions" on "particle" in "position".
 * @tparam field Numeric field for every variable.
 * @param particle Particle number.
 * @param position Particle's position.
 * @param positions Positions of the system's particles.
 * @param r_cut Size of the cut.
 * @return Lennard-Jones potential on "particle".
 */
template<bool tail_correction, typename field>
field LJ_potential(size_t particle, const std::array<field, 3> &position,
                   const Vectors<field> &positions, const SimulationSettings<field> &simulation) {
    field potential{0};
    for (size_t i = 0UL; i < simulation.n_particles; i++) {
        if (particle == i) continue;
        const auto dxij = PBC(positions.e_i[i] - position[0], simulation.box_edge);
        const auto dyij = PBC(positions.e_j[i] - position[1], simulation.box_edge);
        const auto dzij = PBC(positions.e_k[i] - position[2], simulation.box_edge);
        const auto drij2 = dxij * dxij + dyij * dyij + dzij * dzij;
        if (drij2 < simulation.cutoff2) {
            const auto lj1 = field(1) / std::pow(drij2, 6);
            const auto lj2 = field(1) / std::pow(drij2, 3);
            potential += lj1 - lj2;
        }
    }
    potential *= 4;
    if constexpr (tail_correction) return potential + simulation.u_tail_correction;
    else
        return potential;
}
#endif//ESERCIZI_LSN_MS_UTILS_HPP
