//
// Created by Davide Nicoli on 18/07/22.
//

#ifndef ESERCIZI_LSN_MOLECULAR_SYSTEM_HPP
#define ESERCIZI_LSN_MOLECULAR_SYSTEM_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <optional>
#include <random>
#include <type_traits>

#include "algos.hpp"
#include "data_types/measures.hpp"
#include "data_types/settings.hpp"
#include "data_types/vectors.hpp"
#include "estimators/mean.hpp"

using namespace estimators;
namespace fs = std::filesystem;

enum Ensamble { NVT, NVE };

/**
 * Lennard-Jones molecular system
 * @tparam field Numeric field for physical variables
 */
template<typename field, bool tail_corrections, Ensamble ens>
class LJMono {
public:
    typedef field Field;
    LJMono() = delete;
    LJMono(LJMono &) = delete;
    LJMono(const LJMono &) = delete;

    /**
     * Generates a system using settings and positions specified in the respective files.
     * @param settings_path Path to simulation settings.
     * @param positions_path Path to molecular positions.
     */
    LJMono(const fs::path &settings_path, const fs::path &positions_path)
        : m_thermo(settings_path), m_simulation(settings_path),
          m_positions(m_simulation.n_particles, positions_path),
          m_forces(m_simulation.n_particles) {
        m_positions.apply([&](auto &x) { x = pbc(x * m_simulation.box_edge); });
    }

    /**
     * Constructor using settings stored in files. If a velocity file is not provided molecular velocities will be sampled from a Maxwell-Boltzmann distribution.
     * @param settings_path Path to the integrator settings (initial temperature, density, number of particles...)
     * @param positions_path Path to the initial molecular configuration
     * @param velocities_path Optional path to the initial velocity configuration
     * @param rng Random number generator used to eventually sample velocities
     */
    template<class URBG>
    LJMono(const fs::path &settings_path, const fs::path &positions_path,
           const std::optional<fs::path> &velocities_path, URBG &rng)
        : LJMono(settings_path, positions_path) {
        if (velocities_path.has_value()) {
            m_velocities = Vectors<field>(m_simulation.n_particles, velocities_path.value());
            compute_previous_positions(m_positions, m_velocities, m_simulation.delta,
                                       m_simulation.box_edge, m_prev_positions);
        } else {
            m_velocities = Vectors<field>(m_simulation.n_particles);
            generate_velocities(m_velocities, m_thermo.temperature, rng);
            compute_previous_positions(m_positions, m_velocities, m_simulation.delta,
                                       m_simulation.box_edge, m_prev_positions);
        }
    }

    /**
     * Generates a system using settings and positions specified in the respective files. Velocities
     * are generated using a random normal distribution
     * @param settings_path Path to simulation settings.
     * @param positions_path Path to molecular positions.
     * @param rng Random generator.
     */
    template<class URBG, std::enable_if_t<std::is_invocable_v<URBG>, bool> = true>
    LJMono(const fs::path &settings_path, const fs::path &positions_path, URBG &rng)
        : LJMono(settings_path, positions_path) {
        init_velocities(rng);
    }

    /**
     * Generates a system using settings, positions and velocities specified in the respective files
     * @param settings_path Path to simulation settings.
     * @param positions_path Path to molecular positions.
     * @param velocities_path Path to velocities.
     */
    LJMono(const fs::path &settings_path, const fs::path &positions_path,
           const fs::path &velocities_path)
        : LJMono(settings_path, positions_path) {
        init_velocities(velocities_path);
    }


    /**
     * Stores molecules positions in xyz format
     * @param dir Directory where files will be stored
     */
    void save_xyz_positions(const fs::path &dir, size_t frame) const {
        m_positions.save_xyz_configuration(
                (dir / (std::to_string(frame) + ".xyz")).string(),
                [&](const auto x) { return pbc(x) / m_simulation.box_edge; });
    }

    /**
     * Stores molecules positions as a list.
     * @param positions_path Path where positions will be stored.
     */
    void save_positions(const fs::path &positions_path) const {
        m_positions.save_configuration(
                positions_path, [&](const auto x) { return pbc(x) / m_simulation.box_edge; });
    }

    /**
     * Stores molecules positions and velocities as columns of components separated by spaces
     * @param positions_path File where positions will be stored
     * @param velocities_path File where velocities will be stored
     */
    void save_configurations(const fs::path &positions_path,
                             const fs::path &velocities_path) const {
        save_positions(positions_path);
        m_velocities.save_configuration(velocities_path);
    }

    /**
     * Initializes velocities with a normal distribution using the system's thermal properties
     * @param rng Random number generator
     */
    template<class URBG, std::enable_if_t<std::is_invocable_v<URBG>, bool> = true>
    void init_velocities(URBG &rng) {
        m_velocities = Vectors<field>(m_simulation.n_particles);
        generate_velocities(m_velocities, m_thermo.temperature, rng);
        compute_previous_positions(m_positions, m_velocities, m_simulation.delta,
                                   m_simulation.box_edge, m_prev_positions);
    }


    /**
     * Initializes velocities using those stored in the given path.
     * @param velocities_path Path to the velocities.
     */
    void init_velocities(const fs::path &velocities_path) {
        m_velocities = Vectors<field>(m_simulation.n_particles, velocities_path);
        compute_previous_positions(m_positions, m_velocities, m_simulation.delta,
                                   m_simulation.box_edge, m_prev_positions);
    }


    /**
     * Periodic boundary condition using simulator settings
     * @param x A position component
     */
    inline field pbc(field x) const { return PBC(x, m_simulation.box_edge); }


    [[nodiscard]] static constexpr size_t N_VARS() { return 5; }

    /**
     * Names of the thermodynamical variables for which a statistic is computed
     */
    [[nodiscard]] static std::array<std::string, N_VARS()> variable_names() {
        return {"U/N", "K/N", "E/N", "T", "p"};
    }

    [[nodiscard]] size_t time() const { return m_time; }
    void time_step() { m_time++; }

    /**
     * Initializes the radial distribution histogram
     * @param n_bins Number of bins in the histogram.
     */
    void init_radial_func(size_t n_bins) {
        m_n_bins = n_bins;
        m_dr = m_simulation.box_edge / Field(2 * n_bins);
        m_hist.resize(n_bins);
        m_normcoeffs.resize(n_bins);
        m_drs.resize(n_bins);
        const Field k = static_cast<Field>(m_simulation.n_particles) * m_thermo.density;
        Field r{0};
        for (size_t bin = 0; bin < n_bins; bin++) {
            Field dV_r = (4 * M_PI / 3) * (std::pow(r + m_dr, 3) - std::pow(r, 3));
            m_normcoeffs[bin] = Field(1) / (k * dV_r);
            m_drs[bin] = r;
            r += m_dr;
        }
    }

    /**
     * Measures thermodynamical variables at the current time step and computes the forces acting on the molecules
     */
    void measures(Field &potential_energy, Field &kinetic_energy, Field &total_energy,
                  Field &temperature, Field &pressure) {
        field e_pot_{0}, virial_{0};
        std::fill(std::begin(m_forces.e_i), std::end(m_forces.e_i), field(0));
        std::fill(std::begin(m_forces.e_j), std::end(m_forces.e_j), field(0));
        std::fill(std::begin(m_forces.e_k), std::end(m_forces.e_k), field(0));
        for (size_t particle = 0UL; particle + 1 < m_simulation.n_particles; particle++) {
            for (size_t other = particle + 1; other < m_simulation.n_particles; other++) {
                // r_po
                const auto dxij = pbc(m_positions.e_i[particle] - m_positions.e_i[other]);
                const auto dyij = pbc(m_positions.e_j[particle] - m_positions.e_j[other]);
                const auto dzij = pbc(m_positions.e_k[particle] - m_positions.e_k[other]);
                const auto drij2 = dxij * dxij + dyij * dyij + dzij * dzij;
                if (drij2 < m_simulation.cutoff2) {
                    // Computing the addends of both the potential and the virial
                    const auto lj1 = field(1) / std::pow(drij2, 6);
                    const auto lj2 = field(1) / std::pow(drij2, 3);
                    const auto wij = lj1 - lj2 / field(2);
                    // Updating both particles
                    // Force acted on "particle" by "other"
                    const auto fxij = wij * dxij / drij2;
                    const auto fyij = wij * dyij / drij2;
                    const auto fzij = wij * dzij / drij2;
                    m_forces.e_i[particle] += fxij;
                    m_forces.e_j[particle] += fyij;
                    m_forces.e_k[particle] += fzij;
                    // Force acted on "other" by "particle"
                    m_forces.e_i[other] -= fxij;
                    m_forces.e_j[other] -= fyij;
                    m_forces.e_k[other] -= fzij;
                    e_pot_ += (lj1 - lj2);
                    virial_ += wij;
                }
            }
        }
        m_forces *= field(48);
        e_pot_ = 4 * e_pot_ / field(m_simulation.n_particles);
        virial_ *= field(48) / field(3);
        if constexpr (tail_corrections) {
            e_pot_ += m_simulation.u_tail_correction;
            virial_ += m_simulation.W_tail_correction;
        }
        potential_energy = e_pot_;
        if constexpr (ens == NVT) kinetic_energy = field(1.5) * m_thermo.temperature;
        else
            kinetic_energy = m_velocities.full_norm2() / field(2 * m_simulation.n_particles);
        total_energy = potential_energy + kinetic_energy;
        if constexpr (ens == NVT) temperature = m_thermo.temperature;
        else
            temperature = 2 * kinetic_energy / field(3);
        pressure = m_thermo.density * temperature + virial_ / m_simulation.volume;
    }

    /**
     * Measures thermodynamical variables at the current time step and computes the forces acting on the molecules
     * @tparam compute_forces Whether forces should be updated.
     * @tparam vars Variables which will be computed.
     * @param output MeasureOutputs struct to store a compile-time defined number of variables.
     */
    template<bool compute_forces, Variable... vars>
    constexpr void measures(MeasureOutputs<field, vars...> &output) {
        using Outs = MeasureOutputs<field, vars...>;
        constexpr const bool compute_radial = Outs::template has_member<Variable::RadialFn>();
        field e_pot_{0}, virial_{0};
        if constexpr (compute_radial) { std::fill(m_hist.begin(), m_hist.end(), size_t(0)); }
        if constexpr (compute_forces) {
            std::fill(std::begin(m_forces.e_i), std::end(m_forces.e_i), field(0));
            std::fill(std::begin(m_forces.e_j), std::end(m_forces.e_j), field(0));
            std::fill(std::begin(m_forces.e_k), std::end(m_forces.e_k), field(0));
        }
        for (size_t particle = 0UL; particle + 1 < m_simulation.n_particles; particle++) {
            for (size_t other = particle + 1; other < m_simulation.n_particles; other++) {
                // r_po
                const auto dxij = pbc(m_positions.e_i[particle] - m_positions.e_i[other]);
                const auto dyij = pbc(m_positions.e_j[particle] - m_positions.e_j[other]);
                const auto dzij = pbc(m_positions.e_k[particle] - m_positions.e_k[other]);
                const auto drij2 = dxij * dxij + dyij * dyij + dzij * dzij;
                if (drij2 < m_simulation.cutoff2) {
                    // Computing the addends of both the potential and the virial
                    const auto lj1 = field(1) / std::pow(drij2, 6);
                    const auto lj2 = field(1) / std::pow(drij2, 3);
                    const auto wij = lj1 - lj2 / field(2);
                    e_pot_ += (lj1 - lj2);
                    virial_ += wij;
                    if constexpr (compute_forces) {
                        // Updating both particles
                        // Force acted on "particle" by "other"
                        const auto fxij = wij * dxij / drij2;
                        const auto fyij = wij * dyij / drij2;
                        const auto fzij = wij * dzij / drij2;
                        m_forces.e_i[particle] += fxij;
                        m_forces.e_j[particle] += fyij;
                        m_forces.e_k[particle] += fzij;
                        // Force acted on "other" by "particle"
                        m_forces.e_i[other] -= fxij;
                        m_forces.e_j[other] -= fyij;
                        m_forces.e_k[other] -= fzij;
                    }
                }
                if constexpr (compute_radial) {
                    const auto bin = size_t(std::trunc(std::sqrt(drij2) / m_dr));
                    if (bin < m_n_bins) m_hist[bin] += 2;
                }
            }
        }
        if constexpr (compute_forces) m_forces *= field(48);
        e_pot_ = 4 * e_pot_ / field(m_simulation.n_particles);
        virial_ *= (field(48) / field(3));
        if constexpr (tail_corrections) {
            e_pot_ += m_simulation.u_tail_correction;
            virial_ += m_simulation.W_tail_correction;
        }

        output.template push_measure<Variable::PotentialEnergy>(e_pot_);

        field e_kin;
        if constexpr (ens == NVT) e_kin = field(1.5) * m_thermo.temperature;
        else
            e_kin = m_velocities.full_norm2() / field(2 * m_simulation.n_particles);
        output.template push_measure<Variable::KineticEnergy>(e_kin);

        const auto e_tot = e_pot_ + e_kin;
        output.template push_measure<Variable::TotalEnergy>(e_tot);

        field temperature;
        if constexpr (ens == NVT) temperature = m_thermo.temperature;
        else
            temperature = 2 * e_kin / field(3);
        output.template push_measure<Variable::Temperature>(temperature);

        output.template push_measure<Variable::Pressure>(m_thermo.density * temperature +
                                                         virial_ / m_simulation.volume);

        if constexpr (compute_radial) {
            output.template init_vector<Variable::RadialFn>(m_n_bins);
            std::vector<field> g_r(m_hist.size());
            std::transform(m_hist.begin(), m_hist.end(), m_normcoeffs.begin(), g_r.begin(),
                           [](const auto bin, const auto norm) { return field(bin) * norm; });
            output.template push_vector<Variable::RadialFn>(g_r);
        }
    }


    // Thermodynamical configuration: initial temperature and density
    ThermoSettings<field> m_thermo;
    // Simulation settings: number of particles, box edge, periodic volume, ...
    SimulationSettings<field> m_simulation;
    // Molecular positions at the current time step
    Vectors<field> m_positions;
    // Forces acting on each molecule at the current time step
    Vectors<field> m_forces;
    // Molecular velocities at the current time step
    Vectors<field> m_velocities{0};
    // Molecular positions at the previous time step
    Vectors<field> m_prev_positions{0};

    std::vector<field> m_drs{0};

private:
    size_t m_time{0};
    // Constants in the computation of g(r)
    // Number of bins
    size_t m_n_bins{0};
    std::vector<size_t> m_hist{0};
    // dr
    Field m_dr{0};
    // normalization coefficients: 1/(rho*N*dV)
    std::vector<Field> m_normcoeffs{};
};


#endif//ESERCIZI_LSN_MOLECULAR_SYSTEM_HPP
