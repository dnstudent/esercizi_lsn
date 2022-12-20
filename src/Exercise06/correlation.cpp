//
// Created by Davide Nicoli on 05/09/22.
//

/**
 * An executable which computes the auto correlation of one or more series stored in a csv file.
 */

#include <filesystem>

#include <cxxopts.hpp>
#include <rapidcsv.h>

#include "utils.hpp"

namespace co = cxxopts;
namespace fs = std::filesystem;


#define SECTION "06"
#define EXERCISE SECTION "_correlation"

int main(int argc, char const *argv[]) {
    // Defining some flags to control the program output
    cxxopts::Options options("Autocorrelation", "How to run autocorrelation");
    // clang-format off
    options.add_options("Program")
      ("i,in", "Input file", co::value<fs::path>())
      ("o,out", "Output path", co::value<fs::path>())
      ("n,n_lags", "Number of lags to process", co::value<size_t>())
      ("s,skip", "Number of rows to skip from the beginning", co::value<size_t>()->default_value("0"))
      ("h,help", "Print this message")
      ("v,verbose", "Whether to be verbose", co::value<bool>());

    // clang-format on
    auto user_params = options.parse(argc, argv);
    if (user_params.count("help")) {
        std::cout << options.help() << std::endl;
        exit(0);
    }

    const auto INPUT_FILE = user_params["i"].as<fs::path>();
    const auto OUTPUT_FILE = user_params["o"].as<fs::path>();
    // The number of lags for which autocorrelation will be computed.
    const auto N_LAGS = user_params["n"].as<size_t>();
    // The number of samples to skip from the beginning.
    const auto SKIP = user_params["s"].as<size_t>();

    utils::autocorrelation_from<double>(INPUT_FILE, OUTPUT_FILE, N_LAGS, SKIP);

    return 0;
}