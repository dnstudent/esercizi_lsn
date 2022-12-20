#include <array>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <iterator>
#include <valarray>

#include <cxxopts.hpp>
#include <rapidcsv.h>

#include "ariel_random/ariel_random.hpp"
#include "config.hpp"
#include "distributions/cauchy_lorentz.hpp"
#include "estimators/mean.hpp"

#define SECTION "01"
#define EXERCISE SECTION "_3"

using std::string;
using namespace distributions;
using Value = double;
namespace co = cxxopts;
namespace fs = std::filesystem;

typedef std::valarray<Value> point;
typedef std::uniform_real_distribution<Value>::param_type dist_p;

static std::uniform_real_distribution<Value> versor_component_dist(-1, 1);

/**
 * Returns a versor with uniform direction through rejection sampling, to avoid using pi.
 * @param rng A random number generator respecting C++'s UniformRandomBitGenerator
 */
template<class URBG>
Value uniform_versor_y(URBG &rng) {
    point v(2);
    Value x;
    Value y;
    Value norm2;
    do {
        x = versor_component_dist(rng);
        y = versor_component_dist(rng);
        norm2 = x * x + y * y;
    } while (norm2 >= 1);
    return y / std::sqrt(norm2);
}

/**
 * Performs the needle throw
 * @param needle_length Length of the needle.
 * @param a_params A std::uniform_real_distribution<value>::param_type providing the plane where one of the needle ends will be sampled.
 * @param rng A random number generator respecting C++'s UniformRandomBitGenerator.
 * @return The pair of needle extremes.
 */
template<class URBG>
std::pair<Value, Value> throw_needle(Value needle_length,
                                     std::uniform_real_distribution<Value> &y_dist, URBG &rng) {
    const Value a = y_dist(rng);
    const Value delta_y = needle_length * uniform_versor_y(rng);
    const Value b = a + delta_y;
    return {a, b};
}

/**
 * Samples a value for pi simulating Buffon's experiment.
 * @param needle_length The length of the needle.
 * @param half_line_distance Half the distance between two lines.
 * @param n_throws Number of needle throws.
 * @param rng A random number generator respecting C++'s UniformRandomBitGenerator.
 */
template<class URBG>
Value buffon_pi(Value needle_length, Value half_line_distance, size_t n_throws, URBG &rng) {
    size_t n_hits = 0;
    // Letting the first extreme of the needle fall between two lines
    std::uniform_real_distribution<Value> y_dist(-half_line_distance, half_line_distance);
    for (size_t i = 0UL; i < n_throws; i++) {
        auto [a, b] = throw_needle(needle_length, y_dist, rng);
        if ((b > half_line_distance) || (b <= -half_line_distance)) { n_hits++; }
    }
    return needle_length * Value(n_throws) / (Value(n_hits) * half_line_distance);
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
      ("t,n_throws", "Number of needle throws per round", co::value<size_t>()->default_value("1000"))
      ("m,block_size", "Number of blocks for the estimation", co::value<size_t>()->default_value("100"))
      ("n,n_rounds", "Number of rounds", co::value<size_t>()->default_value("100"))
      ("L,needle_length", "Length of the needle", co::value<Value>()->default_value("0.5"))
      ("d,lines_distance", "Distance between the straight lines", co::value<Value>()->default_value("1"));
    // clang-format on
    auto user_params = options.parse(argc, argv);
    if (user_params.count("help")) {
        std::cout << options.help() << std::endl;
        exit(0);
    }

    const Value L = user_params["L"].as<Value>();
    const Value d = user_params["d"].as<Value>();
    const size_t N_THROWS = user_params["t"].as<size_t>();
    const size_t BLOCK_SIZE = user_params["m"].as<size_t>();
    const size_t N_BLOCKS = user_params["n"].as<size_t>();
    const auto OUTPUT_PATH = user_params["o"].as<fs::path>();
    const string PRIMES_SOURCE = user_params["p"].as<string>();
    const size_t PRIMES_LINE = user_params["l"].as<size_t>();
    const string SEEDS_SOURCE = user_params["s"].as<string>();

    if (L >= d) {
        std::cout << "needle_length must be less than lines_distance" << std::endl;
        return 1;
    }

    // ARandom rng(SEEDS_SOURCE, PRIMES_SOURCE, PRIMES_LINE);
    std::minstd_rand rng(1);
    estimators::ProgAvg<Value> mean_estimator;


    const Value d_half = d / Value(2);

    std::vector<Value> pi_samples(BLOCK_SIZE);
    std::vector<Value> mean_estimates(N_BLOCKS);
    auto mit = mean_estimates.begin();
    std::vector<Value> uncert_estimates(N_BLOCKS);
    auto uit = uncert_estimates.begin();
    for (size_t i = 0UL; i < N_BLOCKS; i++) {
        std::generate(pi_samples.begin(), pi_samples.end(),
                      [&]() { return buffon_pi(L, d_half, N_THROWS, rng); });
        std::tie(*mit++, *uit++) = mean_estimator(pi_samples.cbegin(), pi_samples.cend());
    }

    rapidcsv::Document table;
    table.InsertColumn(0, mean_estimates, "estimate");
    table.InsertColumn(1, uncert_estimates, "uncertainty");
    table.RemoveColumn(2);

    if (!fs::exists(OUTPUT_PATH.parent_path())) {
        fs::create_directories(OUTPUT_PATH.parent_path());
    }
    table.Save(OUTPUT_PATH);

    return 0;
}