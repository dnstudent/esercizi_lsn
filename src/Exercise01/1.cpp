//
// Created by Davide Nicoli on 16/03/22.
//

#include <cmath>
#include <filesystem>
#include <functional>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <vector>

#include <cxxopts.hpp>
#include <rapidcsv.h>

#include "ariel_random/ariel_random.hpp"
#include "config.hpp"
#include "distributions/uniform_int.hpp"
#include "estimators/chi2.hpp"
#include "estimators/estimators.hpp"
#include "estimators/mean.hpp"
#include "estimators/variance.hpp"
#include "utils.hpp"

#define SECTION "01"
#define EXERCISE SECTION "_1"

using estimators::ProgAvg;
using estimators::ProgVariance;
using estimators::UniformChi2;
using std::string;
typedef double Value;
typedef std::vector<Value> block_stats;
namespace co = cxxopts;
namespace fs = std::filesystem;

int main(int argc, char *argv[]) {
    // Defining some flags to control the program output
    cxxopts::Options options(EXERCISE, "How to run exercise " EXERCISE);
    // clang-format off
    options.add_options("Program")
      ("out1", "Mean and variance estimates output path", co::value<fs::path>()->default_value(RESULTS_DIR "/" SECTION "/" EXERCISE "_stats.csv"))
      ("out2", "Chi2 statistic output path", co::value<fs::path>()->default_value(RESULTS_DIR "/" SECTION "/"  EXERCISE "_chi.csv"));
    options.add_options("Rng seeding")
      ("p,primes_path", "Prime numbers path", co::value<string>()->default_value(PRIMES_PATH "Primes"))
      ("l,primes_line", "Line in primes_path to use", co::value<size_t>()->default_value("0"))
      ("s, seeds_path", "Seed path", co::value<string>()->default_value(SEEDS_PATH "seed.in"))
      ("h,help", "Print this message");
    options.add_options("Mean and variance")
      ("m,block_size", "Number of samples per block", co::value<size_t>()->default_value("1000"))
      ("n,n_blocks", "Number of blocks", co::value<size_t>()->default_value("100"));
    options.add_options("Pearson X statistic")
      ("n_trials", "Number of trials", co::value<size_t>()->default_value("100"))
      ("n_samples", "Number of samples", co::value<size_t>()->default_value("10000"))
      ("n_intervals", "Number of intervals for the Pearson X statistics", co::value<size_t>()->default_value("100"));
    // clang-format on
    auto user_params = options.parse(argc, argv);
    if (user_params.count("help")) {
        std::cout << options.help() << std::endl;
        exit(0);
    }

    const size_t BLOCK_SIZE = user_params["m"].as<size_t>();
    const size_t N_BLOCKS = user_params["n"].as<size_t>();
    const size_t N_TRIALS = user_params["n_trials"].as<size_t>();
    const size_t N_X_SAMPLES = user_params["n_samples"].as<size_t>();
    const size_t N_INTERVALS = user_params["n_intervals"].as<size_t>();
    const auto EST_OUTPUT_PATH = user_params["out1"].as<fs::path>();
    const auto CHI_OUTPUT_PATH = user_params["out2"].as<fs::path>();
    const string PRIMES_SOURCE = user_params["p"].as<string>();
    const size_t PRIMES_LINE = user_params["l"].as<size_t>();
    const string SEEDS_SOURCE = user_params["s"].as<string>();

    // Building a generator for random numbers between 0 and 1
    ARandom rng(SEEDS_SOURCE, PRIMES_SOURCE, PRIMES_LINE);

    /*----------------------
     * RNG mean and variance
     ----------------------*/
    std::uniform_real_distribution<Value> sampler(0, 1);
    // Buffer for random values
    std::vector<Value> block(BLOCK_SIZE);

    // Vectors storing block statistics
    block_stats mean_estimate(N_BLOCKS), mean_uncert(N_BLOCKS), var_estimate(N_BLOCKS),
            var_uncert(N_BLOCKS);

    // Stateful estimators for block statistics
    ProgAvg<Value> mean_estimator;
    ProgVariance<Value> variance_estimator;
    std::tuple estimators{mean_estimator, variance_estimator};


    std::vector<size_t> blocks(N_BLOCKS);
    std::iota(blocks.begin(), blocks.end(), size_t(0));

    for (const auto block_i: blocks) {
        // Generating a block
        std::generate(block.begin(), block.end(), [&]() { return sampler(rng); });
        // Computing and storing block statistics
        estimators::store_estimates(block.cbegin(), block.cend(), estimators,
                                    utils::snext(mean_estimate.begin(), block_i),
                                    utils::snext(mean_uncert.begin(), block_i),
                                    utils::snext(var_estimate.begin(), block_i),
                                    utils::snext(var_uncert.begin(), block_i));
    }


    // Storing results in a csv document
    rapidcsv::Document table;
    table.SetColumn(0, mean_estimate);
    table.SetColumnName(0, "mean_estimate");
    table.InsertColumn(1, mean_uncert, "mean_error");
    table.InsertColumn(2, var_estimate, "variance_estimate");
    table.InsertColumn(3, var_uncert, "variance_error");

    if (!fs::exists(EST_OUTPUT_PATH.parent_path())) {
        fs::create_directories(EST_OUTPUT_PATH.parent_path());
    }
    table.Save(EST_OUTPUT_PATH);

    /*-----------------------
     * X^2 Pearson statistics
     -----------------------*/
    // vector storing the results
    std::vector<Value> chi2(N_TRIALS);
    // buffer storing each trial's empirical distribution
    std::vector<size_t> O_histogram(N_INTERVALS);

    // expected value for the X^2 statistic
    const Value expected_value = Value(N_X_SAMPLES) / Value(N_INTERVALS);

    distributions::uniform_int<size_t> int_sampler(0, N_INTERVALS - 1);
    std::generate(chi2.begin(), chi2.end(), [&]() {
        std::fill(O_histogram.begin(), O_histogram.end(), 0UL);
        for (size_t i = 0UL; i < N_X_SAMPLES; i++) { O_histogram[int_sampler(rng)]++; }
        return UniformChi2<Value>(expected_value)(O_histogram.cbegin(), O_histogram.cend());
    });

    rapidcsv::Document chi_table;
    chi_table.SetColumn(0, chi2);
    chi_table.SetColumnName(0, "X^2");

    if (!fs::exists(CHI_OUTPUT_PATH.parent_path())) {
        fs::create_directories(CHI_OUTPUT_PATH.parent_path());
    }
    chi_table.Save(CHI_OUTPUT_PATH);

    return 0;
}