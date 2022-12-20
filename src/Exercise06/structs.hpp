//
// Created by Davide Nicoli on 26/08/22.
//

#ifndef ESERCIZI_LSN_06_STRUCTS_HPP
#define ESERCIZI_LSN_06_STRUCTS_HPP

#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <valarray>
#include <vector>

#include <cxxopts.hpp>

namespace fs = std::filesystem;

/**
 * Struct storing program options
 * @tparam var_space The numeric field to use for the variables
 */
template<typename var_space>
struct ExOptions {
    explicit ExOptions(cxxopts::ParseResult &pr)
        : output_dir(pr["out"].as<fs::path>()), n_steps(pr["M"].as<size_t>()),
          block_size(pr["S"].as<size_t>()), warmup_steps(pr["w"].as<size_t>()),
          n_spins(pr["n_spins"].as<size_t>()), metropolis(pr["metropolis"].as<bool>()),
          gibbs(pr["gibbs"].as<bool>()), save_spins(pr["save_spins"].as<bool>()),
          resume(pr["resume"].as<bool>()), J(pr["J"].as<var_space>()), h(pr["B"].as<var_space>()),
          T(pr["T"].as<var_space>()) {
        assert(n_steps % block_size == 0);
    }
    fs::path output_dir;
    size_t n_steps, block_size, n_blocks{n_steps / block_size}, warmup_steps, n_spins;
    bool metropolis, gibbs, save_spins, resume;
    var_space J, h, T;
};

#endif//ESERCIZI_LSN_06_STRUCTS_HPP
