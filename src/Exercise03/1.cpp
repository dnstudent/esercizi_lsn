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


#define SECTION "03"
#define EXERCISE SECTION "_1"

using std::string;
using Value = double;
namespace co = cxxopts;
namespace fs = std::filesystem;

/**
 * Asset price sampler using geometric brownian motion
 */
class GBMAssetSampler {
public:
    /**
     * Instantiate the sampler
     * @param S_0 Initial asset price
     * @param T Maturity
     * @param r Risk-free interest rate
     * @param sigma Price volatility
     */
    GBMAssetSampler(Value S_0, Value T, Value r, Value sigma)
        : m_S_0{S_0}, m_T{T}, m_r{r}, m_sigma{sigma}, m_add1{(r - sigma * sigma / 2) * T},
          m_gauss{0, sigma * std::sqrt(T)} {}


    /**
     * Direct sampler
     * @param rng The random engine used for sampling
     * @return The asset price at time T
     */
    template<class URBG>
    inline Value operator()(URBG &rng) {
        return m_S_0 * std::exp(m_add1 + m_gauss(rng));
    }

    /**
     * Discrete steps sampler
     * @param n_intervals Number of intervals in [0, T]
     * @param rng The random engine used for sampling
     * @return The asset price at time T
     */
    template<class URBG>
    Value operator()(size_t n_intervals, URBG &rng) {
        Value dt = m_T / Value(n_intervals);
        std::normal_distribution<Value> gauss(0, m_sigma * std::sqrt(dt));
        Value add2 = (m_r - m_sigma * m_sigma / 2) * dt;
        Value S_i = m_S_0;
        for (size_t i = 1UL; i <= n_intervals; i++) { S_i *= std::exp(add2 + gauss(rng)); }
        return S_i;
    }

private:
    const Value m_S_0, m_T, m_r, m_sigma, m_add1;
    std::normal_distribution<Value> m_gauss;
};

/**
 * European call option calculator
 */
class CallOptionCalc {
public:
    CallOptionCalc(Value r, Value T, Value K) : m_discount{std::exp(-r * T)}, m_K{K} {}

    /**
     * Computes the european call option price given an underlying asset price at maturity.
     * @param S_T The underlying asset price at maturity
     * @return
     */
    inline Value operator()(Value S_T) const { return m_discount * std::max(S_T - m_K, Value(0)); }

private:
    Value m_discount, m_K;
};

/**
 * European put option calculator
 */
class PutOptionCalc {
public:
    PutOptionCalc(Value r, Value T, Value K) : m_discount{std::exp(-r * T)}, m_K{K} {}

    /**
     * Computes the european put option price given an underlying asset price at maturity.
     * @param S_T The underlying asset price at maturity
     * @return
     */
    inline Value operator()(Value S_T) const { return m_discount * std::max(m_K - S_T, Value(0)); }

private:
    Value m_discount, m_K;
};

/**
 * Utility for block estimation of an option price
 * @param block_size Block size in the estimation procedure
 * @param n_blocks Number of blocks in the estimation procedure
 * @param option_price_calc Option price calculator. Must take an asset price as argument
 * @param asset_sampler Asset sampler. Must take an rng as argument
 * @param table A rapidcsv table in which results will be stored
 * @param section The operation name
 * @param rng The random engine used for sampling
 */
template<typename OptionCost, typename AssetSampler, class URBG>
void estimate_and_store(size_t block_size, size_t n_blocks, OptionCost option_price_calc,
                        AssetSampler asset_sampler, rapidcsv::Document &table,
                        const std::string_view &section, URBG &rng) {
    std::vector<Value> block(block_size);
    estimators::ProgAvg<Value> mean_estimator;
    std::vector<Value> estimates(n_blocks);
    std::vector<Value> uncerts(n_blocks);
    for (size_t i = 0UL; i < n_blocks; i++) {
        std::generate(block.begin(), block.end(),
                      [&]() { return option_price_calc(asset_sampler(rng)); });
        std::tie(estimates[i], uncerts[i]) = mean_estimator(block.cbegin(), block.cend());
    }
    table.InsertColumn(table.GetColumnCount(), estimates, string(section) + "_mean");
    table.InsertColumn(table.GetColumnCount(), uncerts, string(section) + "_error");
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
    options.add_options("Stats")
      ("N,n_blocks", "Number of blocks", co::value<size_t>()->default_value("100"))
      ("M,n_simulations", "Number of simulations", co::value<size_t>()->default_value("10000"))
      ("W,n_intervals", "Number of intervals in [0, T]", co::value<size_t>()->default_value("100"));
    options.add_options("Assets")
      ("S_0", "Initial asset value", co::value<Value>()->default_value("100"))
      ("T,maturity", "Time of expiration", co::value<Value>()->default_value("1"))
      ("r,interest_rate", "Risk-free interest rate", co::value<Value>()->default_value("0.1"))
      ("v,volatility", "Asset price volatility", co::value<Value>()->default_value("0.25"))
      ("K,strike_price", "Strike price", co::value<Value>()->default_value("100"));
    // clang-format on
    auto user_params = options.parse(argc, argv);
    if (user_params.count("help")) {
        std::cout << options.help() << std::endl;
        exit(0);
    }

    const size_t N_SIMULATIONS = user_params["M"].as<size_t>();
    const size_t N_BLOCKS = user_params["N"].as<size_t>();
    const size_t N_INTERVALS = user_params["W"].as<size_t>();
    const auto INITIAL_VALUE = user_params["S_0"].as<Value>();
    const auto MATURITY = user_params["T"].as<Value>();
    const auto INTEREST_RATE = user_params["r"].as<Value>();
    const auto VOLATILITY = user_params["v"].as<Value>();
    const auto STRIKE_PRICE = user_params["K"].as<Value>();
    const fs::path OUTPUT_PATH = user_params["o"].as<fs::path>();
    const string PRIMES_SOURCE = user_params["p"].as<string>();
    const size_t PRIMES_LINE = user_params["l"].as<size_t>();
    const string SEEDS_SOURCE = user_params["s"].as<string>();

    if (N_SIMULATIONS % N_BLOCKS != 0) {
        std::cout
                << "Must choose transform number of blocks which divides the number of simulations"
                << std::endl;
        return 1;
    }
    const size_t BLOCK_SIZE = N_SIMULATIONS / N_BLOCKS;

    rapidcsv::Document table;

    ARandom rng(SEEDS_SOURCE, PRIMES_SOURCE, PRIMES_LINE);

    GBMAssetSampler asset_sampler(INITIAL_VALUE, MATURITY, INTEREST_RATE, VOLATILITY);
    CallOptionCalc call_cost(INTEREST_RATE, MATURITY, STRIKE_PRICE);
    PutOptionCalc put_cost(INTEREST_RATE, MATURITY, STRIKE_PRICE);
    estimate_and_store(BLOCK_SIZE, N_BLOCKS, call_cost, asset_sampler, table, "direct_call", rng);
    estimate_and_store(
            BLOCK_SIZE, N_BLOCKS, call_cost, [&](auto &g) { return asset_sampler(N_INTERVALS, g); },
            table, "discrete_call", rng);
    estimate_and_store(BLOCK_SIZE, N_BLOCKS, put_cost, asset_sampler, table, "direct_put", rng);
    estimate_and_store(
            BLOCK_SIZE, N_BLOCKS, put_cost, [&](auto &g) { return asset_sampler(N_INTERVALS, g); },
            table, "discrete_put", rng);

//    table.RemoveColumn(table.GetColumnCount() - 1);
    if (!fs::exists(OUTPUT_PATH.parent_path())) {
        fs::create_directories(OUTPUT_PATH.parent_path());
    }
    table.Save(OUTPUT_PATH);


    return 0;
}