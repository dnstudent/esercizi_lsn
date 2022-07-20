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
#include "md_integrator/md_integrator.hpp"


#define SECTION "04"
#define EXERCISE SECTION "_1"

using std::string;
using value = double;
namespace co = cxxopts;
namespace fs = std::filesystem;

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
    options.add_options("Stats")
      ("N,n_blocks", "Number of blocks", co::value<size_t>()->default_value("100"))
      ("M,n_simulations", "Number of simulations", co::value<size_t>()->default_value("10000"))
      ("W,n_intervals", "Number of intervals in [0, T]", co::value<size_t>()->default_value("100"));
    options.add_options("Assets")
      ("S_0", "Initial asset value", co::value<value>()->default_value("100"))
      ("T,maturity", "Time of expiration", co::value<value>()->default_value("1"))
      ("r,interest_rate", "Risk-free interest rate", co::value<value>()->default_value("0.1"))
      ("v,volatility", "Asset price volatility", co::value<value>()->default_value("0.25"))
      ("K,strike_price", "Strike price", co::value<value>()->default_value("100"));
    // clang-format on
    auto user_params = options.parse(argc, argv);
    if (user_params.count("help")) {
        std::cout << options.help() << std::endl;
        exit(0);
    }

    const size_t N_SIMULATIONS = user_params["M"].as<size_t>();
    const size_t N_BLOCKS = user_params["N"].as<size_t>();
    const size_t N_INTERVALS = user_params["W"].as<size_t>();
    const auto INITIAL_VALUE = user_params["S_0"].as<value>();
    const auto MATURITY = user_params["T"].as<value>();
    const auto INTEREST_RATE = user_params["r"].as<value>();
    const auto VOLATILITY = user_params["v"].as<value>();
    const auto STRIKE_PRICE = user_params["K"].as<value>();
    const fs::path OUTPUT_PATH = user_params["o"].as<fs::path>();
    const string PRIMES_SOURCE = user_params["p"].as<string>();
    const size_t PRIMES_LINE = user_params["l"].as<size_t>();
    const string SEEDS_SOURCE = user_params["s"].as<string>();

    if (N_SIMULATIONS % N_BLOCKS != 0) {
        std::cout << "Must choose a number of blocks which divides the number of simulations"
                  << std::endl;
        return 1;
    }
    const size_t BLOCK_SIZE = N_SIMULATIONS / N_BLOCKS;

    rapidcsv::Document table;

    ARandom rng(SEEDS_SOURCE, PRIMES_SOURCE, PRIMES_LINE);
    MDIntegrator<value> integrator("", "", rng);
    integrator.thermo_measures();

    return 0;
}