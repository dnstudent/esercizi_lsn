#include <cmath>
#include <filesystem>
#include <iostream>
#include <iterator>
#include <valarray>

#include <cxxopts.hpp>
#include <rapidcsv.h>

#include "ariel_random/ariel_random.hpp"
#include "config.hpp"
#include "distributions/discrete.hpp"
#include "distributions/uniform_angle.hpp"
#include "estimators/mean.hpp"
#include "walkers/walker.hpp"


#define SECTION "02"
#define EXERCISE SECTION "_2"

using std::string;
using Value = double;
namespace co = cxxopts;
namespace fs = std::filesystem;

using discrete_field = int;
using discrete_point = std::valarray<discrete_field>;
using discrete_step = std::valarray<discrete_field>;

using continuum_field = double;
using continuum_point = std::valarray<continuum_field>;
using versor = std::valarray<continuum_field>;

#define ZERO3D                                                                                     \
    { 0, 0, 0 }

/**
 * Factory method for transform walker performing length 1 steps along 3D axes.
 */
auto axial_walker_factory() {
    return Walker{discrete_point(ZERO3D),
                  distributions::UniformDiscrete<discrete_step>(
                          {{0, 0, 1}, {0, 1, 0}, {1, 0, 0}, {0, 0, -1}, {0, -1, 0}, {-1, 0, 0}})};
}

/**
 * Factory method for transform walker performing length 1 steps uniformly in any direcion.
 */
auto omnidir_walker_factory() {
    return Walker{continuum_point(ZERO3D), distributions::Uniform3DDirection<versor>()};
}

/**
 * Utility method to store walk results
 * @param walker_factory The walker factory method to build the walker.
 * @param rng The random number generator to use.
 * @param n_blocks Number of block for block statistics.
 * @param block_size The size of the block for block statistics.
 * @param walk_length The random walk lenght.
 * @param table The rapidcsv::Document to store results to.
 * @param section Name of the operation.
 */
template<typename WalkerFactory, class URBG>
void fill_walk_statistics(WalkerFactory walker_factory, URBG &rng, size_t n_blocks,
                          size_t block_size, size_t walk_length, rapidcsv::Document &table,
                          string &&section) {
    auto walker = walker_factory();
    std::vector<double> distances(block_size);
    estimators::ProgAvg<double> mean_estimator;

    std::vector<double> estimates(n_blocks);
    std::vector<double> uncerts(n_blocks);

    for (size_t i = 0UL; i < n_blocks; i++) {
        std::generate(distances.begin(), distances.end(),
                      [&]() { return utils::norm2(walker.walk(walk_length, rng)); });
        auto [estimate, uncert] = mean_estimator(distances.cbegin(), distances.cend());
        estimates[i] = std::sqrt(estimate);
        uncerts[i] = uncert / (2 * estimates[i]);
        walker.set_current(ZERO3D);
    }

    utils::AppendColumns(table, {section + "_mean", section + "_error"},
                         std::make_tuple(estimates, uncerts));
}

int main(int argc, char const *argv[]) {
    // Defining some flags to control the program output
    cxxopts::Options options(EXERCISE, "How to run exercise " EXERCISE);
    // clang-format off
    options.add_options("Program")
      ("o,out", "Output path", co::value<fs::path>())
      ("h,help", "Print this message");
    options.add_options("Rng seeding")
      ("p,primes_path", "Prime numbers path", co::value<string>()->default_value(PRIMES_PATH "Primes"))
      ("l,primes_line", "Line in primes_path to use", co::value<size_t>()->default_value("1"))
      ("s, seeds_path", "Seed path", co::value<string>()->default_value(SEEDS_PATH "seed.in"));
    options.add_options("Needle")
      ("N,n_blocks", "Number of blocks", co::value<size_t>()->default_value("100"))
      ("M,n_simulations", "Number of simulations", co::value<size_t>()->default_value("10000"))
      ("W,walk_length", "Length of the random walks", co::value<size_t>()->default_value("100"));
    // clang-format on
    auto user_params = options.parse(argc, argv);
    if (user_params.count("help")) {
        std::cout << options.help() << std::endl;
        exit(0);
    }

    const size_t N_SIMULATIONS = user_params["M"].as<size_t>();
    const size_t N_BLOCKS = user_params["N"].as<size_t>();
    const size_t WALK_LENGTH = user_params["W"].as<size_t>();
    const fs::path OUTPUT_PATH = user_params["o"].as<fs::path>();
    const string PRIMES_SOURCE = user_params["p"].as<string>();
    const size_t PRIMES_LINE = user_params["l"].as<size_t>();
    const string SEEDS_SOURCE = user_params["s"].as<string>();

    if (N_SIMULATIONS % N_BLOCKS != 0) {
        std::cout << "Must choose transform number of blocks which divides the number of throws"
                  << std::endl;
        return 1;
    }
    const size_t BLOCK_SIZE = N_SIMULATIONS / N_BLOCKS;

    rapidcsv::Document table;

    ARandom rng(SEEDS_SOURCE, PRIMES_SOURCE, PRIMES_LINE);
    // Producing and storing distance statistics for the axial walker
    fill_walk_statistics(axial_walker_factory, rng, N_BLOCKS, BLOCK_SIZE, WALK_LENGTH, table,
                         "discrete");
    // Producing and storing distance statistics for the uniform direction walker
    fill_walk_statistics(omnidir_walker_factory, rng, N_BLOCKS, BLOCK_SIZE, WALK_LENGTH, table,
                         "continuum");

    if (!fs::exists(OUTPUT_PATH.parent_path())) {
        fs::create_directories(OUTPUT_PATH.parent_path());
    }
    table.Save(OUTPUT_PATH);

    return 0;
}