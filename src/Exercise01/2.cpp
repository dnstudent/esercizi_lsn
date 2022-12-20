#include <array>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <iterator>
#include <valarray>
#include <vector>

#include <cxxopts.hpp>
#include <rapidcsv.h>

#include "ariel_random/ariel_random.hpp"
#include "config.hpp"
#include "distributions/cauchy_lorentz.hpp"
#include "distributions/exponential.hpp"
#include "distributions/uniform_int.hpp"

#define SECTION "01"
#define EXERCISE SECTION "_2"

using std::string;
using namespace distributions;
typedef double Value;
namespace co = cxxopts;
namespace fs = std::filesystem;


template<typename It, class Distribution, class URBG>
void sample_averages(It first, It last, size_t N, Distribution &dist, URBG &rng) {
    std::valarray<typename Distribution::result_type> buffer(N);
    std::generate(first, last, [&]() {
        std::generate(std::begin(buffer), std::end(buffer), [&]() { return dist(rng); });
        return Value(buffer.sum()) / Value(buffer.size());
    });
}


int main(int argc, char *argv[]) {
    // Defining some flags to control the program output
    cxxopts::Options options(EXERCISE, "How to run exercise " EXERCISE);
    // clang-format off
    options.add_options("Program")
      ("o,out", "Output path", co::value<fs::path>()->default_value(RESULTS_DIR "/" SECTION "/" EXERCISE ".csv"));
    options.add_options("Rng seeding")
      ("p,primes_path", "Prime numbers path", co::value<string>()->default_value(PRIMES_PATH "Primes"))
      ("l,primes_line", "Line in primes_path to use", co::value<size_t>()->default_value("0"))
      ("s, seeds_path", "Seed path", co::value<string>()->default_value(SEEDS_PATH "seed.in"))
      ("h,help", "Print this message");
    options.add_options("Dices")
      ("n,n_realizations", "Number of realizations of the dices throws", co::value<size_t>()->default_value("10000"));
    // clang-format on
    auto user_params = options.parse(argc, argv);
    if (user_params.count("help")) {
        std::cout << options.help() << std::endl;
        exit(0);
    }

    const size_t N_REALIZATIONS = user_params["n"].as<size_t>();
    const auto OUTPUT_PATH = user_params["o"].as<fs::path>();
    const string PRIMES_SOURCE = user_params["p"].as<string>();
    const size_t PRIMES_LINE = user_params["l"].as<size_t>();
    const string SEEDS_SOURCE = user_params["s"].as<string>();

    ARandom rng(SEEDS_SOURCE, PRIMES_SOURCE, PRIMES_LINE);
    // Declaring the dices
    uniform_int<unsigned short> uniform_dice(1, 6);
    exponential<Value> exponential_dice(1.);
    cauchy_lorentz<Value> lorentzian_dice(0, 1);
    std::array<size_t, 4> Ns{1, 2, 10, 100};

    rapidcsv::Document table;
    size_t col_index = 1;
    for (const auto N: Ns) {
        std::vector<Value> realizations(N_REALIZATIONS);
        sample_averages(realizations.begin(), realizations.end(), N, uniform_dice, rng);
        table.InsertColumn(col_index, realizations, "uniform_" + std::to_string(N));
        col_index++;
        sample_averages(realizations.begin(), realizations.end(), N, exponential_dice, rng);
        table.InsertColumn(col_index, realizations, "exponential_" + std::to_string(N));
        col_index++;
        sample_averages(realizations.begin(), realizations.end(), N, lorentzian_dice, rng);
        table.InsertColumn(col_index, realizations, "lorentzian_" + std::to_string(N));
        col_index++;
    }
    table.RemoveColumn(0);
    if (!fs::exists(OUTPUT_PATH.parent_path())) {
        fs::create_directories(OUTPUT_PATH.parent_path());
    }
    table.Save(OUTPUT_PATH);
    return 0;
}