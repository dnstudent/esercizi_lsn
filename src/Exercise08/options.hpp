//
// Created by Davide Nicoli on 25/10/22.
//

#ifndef ESERCIZI_LSN_08_OPTIONS_HPP
#define ESERCIZI_LSN_08_OPTIONS_HPP

#include <filesystem>

#include <cxxopts.hpp>

namespace fs = std::filesystem;

namespace ex08 {
    template<typename field>
    struct ExOptions {
        explicit ExOptions(cxxopts::ParseResult &pr)
            : out(pr["out"].as<fs::path>()), n_T_steps(pr["N"].as<size_t>()),
              n_explore_steps(pr["W"].as<size_t>()), n_blocks(pr["n_blocks"].as<size_t>()),
              block_size(pr["block_size"].as<size_t>()), T0(pr["T0"].as<field>()),
              Tf(pr["Tf"].as<field>()), stddev(pr["stddev"].as<field>()) {
            if (!fs::exists(out)) fs::create_directories(out);
            const auto params = pr["p0"].as<std::vector<field>>();
            if (params.size() != 2)
                throw std::runtime_error(
                        "Wrong number of parameters in the initial guess. Must be 2.");
            m0 = params[0];
            s0 = params[1];
        };

        fs::path out;
        size_t n_T_steps, n_explore_steps, n_blocks, block_size;
        field T0, Tf, stddev, m0{}, s0{};
    };

    template<typename field>
    struct ExPsiOptions {
        explicit ExPsiOptions(cxxopts::ParseResult &pr)
            : out(pr["out"].as<fs::path>()), n_blocks(pr["n_blocks"].as<size_t>()),
              block_size(pr["block_size"].as<size_t>()), n_bins(pr["n_bins"].as<size_t>()),
              mu(pr["mu"].as<field>()), sigma(pr["sigma"].as<field>()) {
            if (!fs::exists(out)) fs::create_directories(out);
            const auto bounds_in = pr["bounds"].as<std::vector<field>>();
            if (bounds_in.size() != 2)
                throw std::runtime_error("'bounds' must be a length 2 vector: --bounds=a,b");
            a = bounds_in[0];
            b = bounds_in[1];
        };

        fs::path out;
        size_t n_blocks, block_size, n_bins;
        field mu, sigma, a{}, b{};
    };
}// namespace ex08

#endif//ESERCIZI_LSN_08_OPTIONS_HPP
