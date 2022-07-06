//
// Created by Davide Nicoli on 16/03/22.
//

#include <algorithm>
#include <cmath>
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
#include "estimators/estimators.hpp"
#include "estimators/mean.hpp"
#include "estimators/variance.hpp"
#include "utils.hpp"

//#define BLOCK_SIZE 1000

using estimators::StatefulMean;
using estimators::StatefulVariance;
using std::string;
typedef double value;
typedef std::vector<value> block_stats;
namespace co = cxxopts;

int main(int argc, char *argv[]) {
    // Defining some flags to control the program output
    cxxopts::Options options("01_1", "How to run exercise 01.1");
    // clang-format off
    options.add_options("Program")
      ("o,out", "Output path", co::value<string>()->default_value(RESULTS_DIR "01_1.csv"));
    options.add_options("Rng seeding")
      ("p,primes_path", "Prime numbers path", co::value<string>()->default_value(PRIMES_PATH "Primes"))
      ("l,primes_line", "Line in primes_path to use", co::value<size_t>()->default_value("0"))
      ("s, seeds_path", "Seed path", co::value<string>()->default_value(SEEDS_PATH "seed.in"))
      ("h,help", "Print this message");
    options.add_options("Mean and variance")
      ("m,block_size", "Number of samples per block", co::value<size_t>()->default_value("1000"))
      ("n,n_blocks", "Number of blocks", co::value<size_t>()->default_value("100"));
    options.add_options("Pearson X statistic")
      ("N,n_samples", "Number of samples", co::value<size_t>()->default_value("1000000"))
      ("k,n_intervals", "Number of intervals for the Pearson X statistics", co::value<size_t>()->default_value("100"));
    // clang-format on
    auto user_params = options.parse(argc, argv);
    if (user_params.count("help")) {
        std::cout << options.help() << std::endl;
        exit(0);
    }

    const size_t BLOCK_SIZE = user_params["m"].as<size_t>();
    const size_t N_BLOCKS = user_params["n"].as<size_t>();
    const size_t N_X_SAMPLES = user_params["N"].as<size_t>();
    const size_t N_INTERVALS = user_params["k"].as<size_t>();
    const string OUTPUT_PATH = user_params["o"].as<string>();
    const string PRIMES_SOURCE = user_params["p"].as<string>();
    const size_t PRIMES_LINE = user_params["l"].as<size_t>();
    const string SEEDS_SOURCE = user_params["s"].as<string>();

    // Building a generator for random numbers between 0 and 1
    ARandom rng(SEEDS_SOURCE, PRIMES_SOURCE, PRIMES_LINE);

    /*----------------------
     * RNG mean and variance
     ----------------------*/
    std::uniform_real_distribution<value> sampler(0, 1);
    // Buffer for random values
    std::vector<value> block(BLOCK_SIZE);

    // Vectors storing block statistics
    block_stats mean_estimate(N_BLOCKS), mean_uncert(N_BLOCKS), var_estimate(N_BLOCKS),
            var_uncert(N_BLOCKS);

    // Stateful estimators for block statistics
    StatefulMean<value> mean_estimator;
    StatefulVariance<value> variance_estimator(0.5);
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
    table.SetColumn(0, blocks);
    table.SetColumnName(0, "block");
    table.InsertColumn(1, mean_estimate, "mean_estimate");
    table.InsertColumn(2, mean_uncert, "mean_error");
    table.InsertColumn(3, var_estimate, "variance_estimate");
    table.InsertColumn(4, var_uncert, "variance_error");
    table.Save(OUTPUT_PATH);

    /*-----------------------
     * X^2 Pearson statistics
     -----------------------*/
    std::mt19937_64 g;
    std::uniform_int_distribution<size_t> int_sampler(0, N_INTERVALS - 1);
    std::vector<size_t> O_histogram(N_INTERVALS, 0);
    for (size_t i = 0UL; i < N_X_SAMPLES; i++) {
        //        O_histogram[size_t(std::floor(rng.Rannyu(0, double(N_INTERVALS))))]++;
        O_histogram[int_sampler(rng)]++;
    }
    const value expected_value = value(N_X_SAMPLES) / value(N_INTERVALS);
    //    const value variance =
    //            value(N_X_SAMPLES) * (value(N_INTERVALS - 1) / value(N_INTERVALS)) / value(N_INTERVALS);
    auto chi2 = std::transform_reduce(O_histogram.cbegin(), O_histogram.cend(), value(0),
                                      std::plus<>(), [=](const auto Oi) {
                                          auto num = value(Oi) - expected_value;
                                          return num * num / expected_value;
                                      });

    std::cout << "Pearson's X statistic was " << chi2 << " against a number of intervals of "
              << N_INTERVALS << std::endl;

    return 0;
}