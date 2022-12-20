//
// Created by Davide Nicoli on 28/09/22.
//

#ifndef ESERCIZI_LSN_MS_SETTINGS_HPP
#define ESERCIZI_LSN_MS_SETTINGS_HPP

#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>

namespace fs = std::filesystem;

/**
 * Collection of thermodynamic settings
 * @tparam field
 */
template<typename field>
struct ThermoSettings {
    ThermoSettings(field T, field rho) : temperature{T}, density{rho} {}

    explicit ThermoSettings(const fs::path &settings_path) {
        std::ifstream settings(settings_path);
        if (!settings.is_open())
            throw std::runtime_error("Could not open " + settings_path.string());
        bool _bucket, _resume;
        field T, rho;
        size_t _n_part;
        settings >> _bucket >> _resume >> T >> _n_part >> rho;
        settings.close();
        *this = ThermoSettings<field>(T, rho);
    }

    field temperature, density;
};

/**
 * Collection of MD related settings
 * @tparam field
 */
template<typename field>
struct SimulationSettings {
    SimulationSettings(size_t n_part, size_t n_b, size_t n_s, field r, field d, field density/*,
                       bool is_mc*/)
        : n_particles{n_part}, n_blocks{n_b},
          block_size{n_s}, cutoff{r}, cutoff2{r * r}, volume{field(n_particles) / density},
          box_edge{std::cbrt(volume)}, delta{d}, delta2{delta * delta}, dbldelta{2 * delta},
          u_tail_correction{8 * M_PI * density *
                            (1.0 / (9 * std::pow(cutoff, 9)) - 1.0 / (3 * std::pow(cutoff, 3)))},
          W_tail_correction{field(96) * M_PI * density * static_cast<field>(n_particles) *
                            (1.0 / (9 * std::pow(cutoff, 9)) - 1.0 / (6 * std::pow(cutoff, 3)))} {}

    explicit SimulationSettings(const fs::path &settings_path) {
        std::ifstream settings(settings_path);
        if (!settings.is_open())
            throw std::runtime_error("Could not open " + settings_path.string());
        bool is_mc, _resume;
        field T, rho, r, deltat;
        size_t n_part, n_b, n_s;
        settings >> is_mc >> _resume >> T >> n_part >> rho >> r >> deltat >> n_b >> n_s;
        settings.close();
        *this = SimulationSettings<field>(n_part, n_b, n_s, r, deltat, rho /*, is_mc*/);
    }

    size_t n_particles{}, n_blocks{}, block_size{};
    field cutoff, cutoff2, volume, box_edge, delta, delta2, dbldelta;
    // tail corrections;
    field u_tail_correction, W_tail_correction;
};


#endif//ESERCIZI_LSN_MS_SETTINGS_HPP
