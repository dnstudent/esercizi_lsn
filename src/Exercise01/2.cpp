#include <cmath>
#include <iostream>
#include <iterator>
#include <vector>

#include "argparser.hpp"
#include "config.hpp"
#include "distributions/exponential.hpp"
#include "parser_opts.hpp"
#include "random.hpp"
#include "random_wrapper.old.hpp"

using std::next;
using std::string;
using namespace distributions;
using value = double;
using block_stat = std::vector<value>;
namespace ap = stypox;
using ap::ArgParser;
namespace fs = std::filesystem;

int main(int argc, char const *argv[]) {
  bool help = false;
  string output_path = "simulation_01.2.csv";
  string seed_source = SEEDS_PATH "seed.in";
  string primes_source = PRIMES_PATH "Primes";
  size_t n_intervals = 1E2, n_samples = 1E2, sample_size = 1E6;

  ArgParser p{
      std::make_tuple(
          stypox::HelpSection{"\nHelp options:"},
          stypox::SwitchOption{"help", help, stypox::args("-?", "-h", "--help"),
                               "Show help screen"},
          ap::Option{"n_samples", n_samples, ap::args("--n_samples="),
                     "Number of random samples"},
          ap::Option{"sample_size", sample_size, ap::args("--sample_size="),
                     "Number of blocks to split the samples into"},
          ap::Option{"n_intervals", n_intervals, ap::args("--n_intervals="),
                     "Number of cuts in the [0, 1) interval"},
          ap::Option{
              "output", output_path, ap::args("-o=", "--output="),
              "Path where the csv-formatted table of results will be stored"},
          rng_options::primes_path(primes_source),
          rng_options::seed_path(seed_source)),
      "Exercise 01.1.3"};
  p.parse(argc, argv);
  if (help) {
    std::cout << p.help() << std::endl;
    return 0;
  }
  p.validate();

  ARandom rng(seed_source, primes_source);

  exponential exp;

  lsn_rng::close(&rng);
}