#include <cmath>
#include <filesystem>
#include <iostream>
#include <iterator>
#include <random>
#include <valarray>

#include <cxxopts.hpp>
#include <rapidcsv.h>

#include "ariel_random/ariel_random.hpp"
#include "config.hpp"
#include "estimators/mean.hpp"


#define SECTION "02"
#define EXERCISE SECTION "_1"

using std::string;
using Value = double;
namespace co = cxxopts;
namespace fs = std::filesystem;

static const Value M_SQRT_2_PI = std::sqrt(2 * M_PI);

/**
 * The exercise integrand sampling with transform uniform distribution.
 * @param x x
 */
Value f1(const Value x) { return M_PI_2 * cos(x * M_PI_2); }

/**
 * The exercise integrand sampling with transform normal distribution.
 * @param x x
 * @param mu Normal distribution's expected value.
 * @param sigma Normal distribution's standard deviation.
 */
Value f2(const Value x, const Value mu, const Value sigma) {
    // Using a simmetric domain and a /2 factor for better performances, as both the integrand and the gaussian are simmetric.
    if (-1 < x && x < 1) {
        const auto y = (x - mu) / sigma;
        return f1(x) * std::exp(y * y / 2) * M_SQRT_2_PI * sigma / 2;
    } else
        return 0.;
}

int main(int argc, char const *argv[]) {
    // Defining some flags to control the program output
    cxxopts::Options options(EXERCISE, "How to run exercise " EXERCISE);
    // clang-format off
    options.add_options("Program")
      ("o,out", "Output path", co::value<fs::path>()->default_value(RESULTS_DIR "/" SECTION "/" EXERCISE ".csv"))
      ("h,help", "Print this message");
    options.add_options("Rng seeding")
      ("p,primes_path", "Prime numbers path", co::value<string>()->default_value(PRIMES_PATH "Primes"))
      ("l,primes_line", "Line in primes_path to use", co::value<size_t>()->default_value("1"))
      ("s, seeds_path", "Seed path", co::value<string>()->default_value(SEEDS_PATH "seed.in"));
    options.add_options("Needle")
      ("N,n_blocks", "Number of blocks", co::value<size_t>()->default_value("1000"))
      ("M,n_throws", "Number of throws", co::value<size_t>()->default_value("10000"));
    // clang-format on
    auto user_params = options.parse(argc, argv);
    if (user_params.count("help")) {
        std::cout << options.help() << std::endl;
        exit(0);
    }

    const size_t N_THROWS = user_params["M"].as<size_t>();
    const size_t N_BLOCKS = user_params["N"].as<size_t>();
    const fs::path OUTPUT_PATH = user_params["o"].as<fs::path>();
    const string PRIMES_SOURCE = user_params["p"].as<string>();
    const size_t PRIMES_LINE = user_params["l"].as<size_t>();
    const string SEEDS_SOURCE = user_params["s"].as<string>();

    if (N_THROWS % N_BLOCKS != 0) {
        std::cout << "Must choose transform number of blocks which divides the number of throws"
                  << std::endl;
        return 1;
    }
    const size_t BLOCK_SIZE = N_THROWS / N_BLOCKS;

    ARandom rng(SEEDS_SOURCE, PRIMES_SOURCE, PRIMES_LINE);
    rapidcsv::Document table;
    size_t col_index = 0;
    // The block
    std::vector<Value> sample(BLOCK_SIZE);
    // Statistics storages
    std::vector<Value> integral_estimates(N_BLOCKS);
    std::vector<Value> integral_uncert(N_BLOCKS);
    // The block estimator
    estimators::ProgAvg<Value> integral_estimator;

    // .1 Sample using the uniform distribution.
    std::uniform_real_distribution<Value> unit_dist(0, 1);
    for (size_t i = 0UL; i < N_BLOCKS; i++) {
        std::generate(sample.begin(), sample.end(), [&]() { return f1(unit_dist(rng)); });
        std::tie(integral_estimates[i], integral_uncert[i]) =
                integral_estimator(sample.cbegin(), sample.cend());
    }
    table.InsertColumn(col_index++, integral_estimates, "I1_estimate");
    table.InsertColumn(col_index++, integral_uncert, "I1_error");


    // .2
    std::normal_distribution<Value> gauss_dist(0, 1);
    for (size_t i = 0UL; i < N_BLOCKS; i++) {
        std::generate(sample.begin(), sample.end(), [&]() { return f2(gauss_dist(rng), 0, 1); });
        std::tie(integral_estimates[i], integral_uncert[i]) =
                integral_estimator(sample.cbegin(), sample.cend());
    }
    table.InsertColumn(col_index++, integral_estimates, "I2_estimate");
    table.InsertColumn(col_index++, integral_uncert, "I2_error");
    table.RemoveColumn(col_index);

    if (!fs::exists(OUTPUT_PATH.parent_path())) {
        fs::create_directories(OUTPUT_PATH.parent_path());
    }
    table.Save(OUTPUT_PATH);
    return 0;
}