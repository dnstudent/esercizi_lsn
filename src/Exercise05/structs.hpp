//
// Created by Davide Nicoli on 26/08/22.
//

#ifndef ESERCIZI_LSN_05_STRUCTS_HPP
#define ESERCIZI_LSN_05_STRUCTS_HPP

#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <valarray>
#include <vector>

#include <cxxopts.hpp>

namespace fs = std::filesystem;

template<typename real_space>
struct ExOptions {
    explicit ExOptions(cxxopts::ParseResult &pr)
        : output_dir(pr["out"].as<fs::path>()), n_throws(pr["M"].as<size_t>()),
          n_blocks(pr["N"].as<size_t>()), warmup_steps(pr["w"].as<size_t>()),
          sample_uniform(pr["uniform"].as<bool>()), sample_gauss(pr["gauss"].as<bool>()),
          sample_s(pr["orbital_s"].as<bool>()), sample_2p(pr["orbital_2p"].as<bool>()),
          save_positions(pr["positions"].as<bool>()) {
        if (n_throws % n_blocks != 0) {
            throw std::runtime_error(
                    "The number of throws must be divisible by the number of blocks.");
        }
        const auto step_params = pr["steppers_config"].as<std::vector<real_space>>();
        if (step_params.size() != 4) {
            throw cxxopts::option_syntax_exception(
                    "The stepper configuration must be a list of four scalars.");
        }
        step_unif_s = step_params[0];
        step_gauss_s = step_params[1];
        step_unif_2p = step_params[2];
        step_gauss_2p = step_params[3];


        const auto S0_v = pr["starting_point"].as<std::vector<real_space>>();
        S0 = std::valarray<real_space>(S0_v.size());
        std::copy(S0_v.cbegin(), S0_v.cend(), std::begin(S0));
    }
    fs::path output_dir;
    size_t n_throws, n_blocks, block_size{n_throws / n_blocks}, warmup_steps;
    bool sample_uniform, sample_gauss, sample_s, sample_2p, save_positions;
    real_space step_unif_s, step_gauss_s, step_unif_2p, step_gauss_2p;
    std::valarray<real_space> S0;
};

#endif//ESERCIZI_LSN_05_STRUCTS_HPP
