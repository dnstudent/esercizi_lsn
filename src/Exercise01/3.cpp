//
// Created by Davide Nicoli on 16/03/22.
//

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <numeric>
#include <vector>

#include "argparser.hpp"
#include "config.hpp"
#include "parser_opts.hpp"
#include "random.hpp"
#include "random_wrapper.old.hpp"
#include "rapidcsv.h"

using std::next;
using std::string;
using value = double;
using block_stat = std::vector<value>;
namespace ap = stypox;
using ap::ArgParser;
namespace fs = std::filesystem;

template <typename It, typename Value>
value chi_squared(It counts_first, It counts_last, const Value expected) {
  return std::accumulate(counts_first, counts_last, Value(0),
                         [&](const Value acc, const auto oi) {
                           return acc + pow(Value(oi) - expected, 2);
                         }) /
         expected;
}

int main(int argc, char const *argv[]) {
  bool help = false;
  string output_path = "simulation_01.1.3.csv";
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
  lsn_rng::seed(&rng, seed_source, primes_source);

  std::vector<std::vector<size_t>> counts(n_intervals,
                                          std::vector<size_t>(n_samples, 0));
  std::vector<size_t> sample(sample_size);
  for (size_t i = 0; i < n_samples; i++) {
    for (size_t j = 0; j < sample_size; j++) {
      auto val = size_t(floor(rng.Rannyu(0.0, value(n_intervals))));
      counts[val][i]++;
    }
  }
  lsn_rng::close(&rng);

  std::vector<value> chi(n_intervals, 0);
  const value expected = value(sample_size) / value(n_intervals);

  std::transform(counts.cbegin(), counts.cend(), chi.begin(),
                 [&](const std::vector<size_t> &bin_counts) {
                   return chi_squared(bin_counts.cbegin(), bin_counts.cend(),
                                      expected);
                 });

  std::vector<double> x(n_intervals);
  for (size_t i = 0; i < n_intervals; i++)
    x[i] = double(i) / double(n_intervals);

  // Storing results in a csv table
  rapidcsv::Document table;
  table.SetColumn(0, x);
  table.SetColumnName(0, "x");
  table.InsertColumn(1, chi, "chi_squared");

  //    table.InsertColumn(3, try_est, "error2");

  table.Save(output_path);
  return 0;
}