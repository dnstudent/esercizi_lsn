//
// Created by Davide Nicoli on 18/07/22.
//

#ifndef ESERCIZI_LSN_MD_INTEGRATOR_HPP
#define ESERCIZI_LSN_MD_INTEGRATOR_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <optional>
#include <random>
#include <string_view>
#include <tuple>
#include <valarray>

#include "estimators/mean.hpp"
#include "structures.hpp"


using namespace estimators;


template<typename field>
class MDIntegrator {
public:
    MDIntegrator() = delete;
    MDIntegrator(MDIntegrator &) = delete;
    MDIntegrator(const MDIntegrator &) = delete;
    MDIntegrator(std::string_view settings_path, std::string_view configuration_path,
                 std::string_view velocities_path)
        : MDIntegrator(settings_path, configuration_path) {
        m_velocities = Vectors<field>(m_simulation.n_particles, velocities_path);
        set_prev_positions();
    }
    template<typename URBG>
    MDIntegrator(std::string_view settings_path, std::string_view configuration_path, URBG &rng)
        : MDIntegrator(settings_path, configuration_path) {
        m_velocities = Vectors<field>(m_simulation.n_particles);
        generate_velocities(rng);
        set_prev_positions();
    }

    /*    void load_configuration(std::string_view configuration_path) {
            m_positions = Vectors<field>(m_simulation.n_particles, configuration_path);
        }
        void save_configurations(std::string_view positions_path,
                                 std::string_view velocities_path) const {
            m_positions.save_configuration(positions_path, [&](const auto x) { return x; });
            m_velocities.save_configuration(velocities_path);
        }
        void load_velocities(std::string_view velocities_path) {
            m_velocities = Vectors<field>(m_simulation.n_particles, velocities_path);
        }*/

    void step() {

    }

    std::valarray<field> thermo_measures() {
        field e_pot{0}, virial{0};
        for (size_t i = 0UL; i < m_simulation.n_particles - 1; i++) {
            for (size_t j = i + 1; j < m_simulation.n_particles; j++) {
                const auto dx = m_positions.at<0>(i) - m_positions.at<0>(j);
                const auto dy = m_positions.at<1>(i) - m_positions.at<0>(1);
                const auto dz = m_positions.at<2>(i) - m_positions.at<0>(2);
                const auto dr = std::hypot(dx, dy, dz);
                if (dr < m_thermo.cutoff) {
                    auto lj1 = field{1} / std::pow(dr, 12);
                    auto lj2 = field{1} / std::pow(dr, 6);
                    e_pot += lj1 + lj2;
                    virial += lj1 - lj2 / field{2};
                }
            }
        }
        e_pot *= field{4};
        virial *= (field{48} / field{3});
        field e_kin = m_velocities.full_norm2() / field(2);
        field temperature = 2 * e_kin / field(3 * m_simulation.n_particles);
        field pressure = m_thermo.density * temperature + virial / m_simulation.volume;
        return {e_pot, e_kin, e_pot + e_kin, temperature, pressure};
    }

    field pbc(field x) { return x - m_simulation.box_edge * std::rint(x / m_simulation.box_edge); }

private:
    MDIntegrator(std::string_view settings_path, std::string_view configuration_path)
        : m_thermo(settings_path), m_simulation(settings_path),
          m_positions(m_simulation.n_particles, configuration_path) {
        m_positions.transform([&](auto x) { return pbc(x * m_simulation.box_edge); });
    }

    ThermoSettings<field> m_thermo;
    SimulationSettings<field> m_simulation;
    Vectors<field> m_positions;
    Vectors<field> m_velocities{0};
    Vectors<field> m_prev_positions{0};

    template<typename URBG>
    void generate_velocities(URBG &rng) {
        std::normal_distribution<field> gauss(0, std::sqrt(m_thermo.temperature));
        std::generate(std::begin(m_velocities.e_i), std::end(m_velocities.e_i),
                      [&]() { return gauss(rng); });
        std::generate(std::begin(m_velocities.e_j), std::end(m_velocities.e_j),
                      [&]() { return gauss(rng); });
        std::generate(std::begin(m_velocities.e_k), std::end(m_velocities.e_k),
                      [&]() { return gauss(rng); });
        m_velocities -= m_velocities.mean();
        field v_scale_factor =
                std::sqrt(3 * m_thermo.temperature /
                          (m_velocities.full_norm2() / field(m_simulation.n_particles)));
        m_velocities *= v_scale_factor;
    }

    void set_prev_positions() {
        // convoluted approach to compute previous positions, but it avoids multiple allocations
        // without writing complex code
        m_prev_positions = m_velocities;
        m_prev_positions *= m_simulation.dt;
        m_prev_positions -= m_positions;
        m_prev_positions.transform([&](auto x) { return pbc(-x); });
    }

    inline field reduced_lj_potential(field dr) {
        return field(1.0) / std::pow(dr, 12) - field(1.0) / std::pow(dr, 6);
    }

    inline field red_lj_virial(field dr) {}
};

#endif//ESERCIZI_LSN_MD_INTEGRATOR_HPP
