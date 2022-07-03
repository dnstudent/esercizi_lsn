//
// Created by Davide Nicoli on 16/03/22.
//

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>
#include <vector>

#include "argparser.hpp"
#include "config.hpp"
#include "estimators/estimators.hpp"
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

int main(int argc, char const *argv[]) {
  bool help = false;
  string output_path = "simulation_01.1.csv";
  string seed_source = SEEDS_PATH "seed.in";
  string primes_source = PRIMES_PATH "Primes";
  size_t SAMPLE_SIZE = 1E6, N_BLOCKS = 1E2;

  ArgParser p{
      std::make_tuple(
          stypox::HelpSection{"\nHelp options:"},
          stypox::SwitchOption{"help", help, stypox::args("-?", "-h", "--help"),
                               "Show help screen"},
          ap::Option{"n_throws", SAMPLE_SIZE, ap::args("-M=", "--sample_size="),
                     "Number of random samples"},
          ap::Option{"n_blocks", N_BLOCKS, ap::args("-N=", "--n_blocks="),
                     "Number of blocks to split the samples into"},
          ap::Option{
              "output", output_path, ap::args("-o=", "--output="),
              "Path where the csv-formatted table of results will be stored"},
          rng_options::primes_path(primes_source),
          rng_options::seed_path(seed_source)),
      "Exercise 01.1 .1 e .2"};
  p.parse(argc, argv);
  if (help) {
    std::cout << p.help() << std::endl;
    return 0;
  }
  p.validate();
  if (SAMPLE_SIZE % N_BLOCKS != 0) {
    std::cout << "n_throws must be divisible by n_blocks" << std::endl;
    return 1;
  }

  ARandom rng(seed_source, primes_source);

  // Generating random values
  std::vector<value> sample(SAMPLE_SIZE);
  std::generate(sample.begin(), sample.end(),
                [&rng]() { return rng.Rannyu(); });
  lsn_rng::close(&rng);

  const size_t L = SAMPLE_SIZE / N_BLOCKS;

  // 1.
  const estimators::blocks::mean::Arithmetic<value> mean_estimator(L);
  const block_stat m_estimates =
      mean_estimator.estimate(sample.cbegin(), sample.cend());
  const block_stat m_errors =
      mean_estimator.std(sample.cbegin(), sample.cend());

  // 2.
  const estimators::blocks::mean::Arithmetic<value> variance_estimator(L);
  const auto v_estimates = variance_estimator.estimate(
      sample.cbegin(), sample.cend(),
      [](const value r) { return pow(r - 0.5, 2); });
  const auto v_errors =
      variance_estimator.std(sample.cbegin(), sample.cend(),
                             [](const value r) { return pow(r - 0.5, 2); });

  std::vector<size_t> blocks(N_BLOCKS);
  for (size_t i = 0; i < blocks.size(); i++)
    blocks[i] = i * L;

  // Storing results in a csv table
  rapidcsv::Document table;
  table.SetColumn(0, blocks);
  table.SetColumnName(0, "block");
  table.InsertColumn(1, m_estimates, "mean_estimate");
  table.InsertColumn(2, m_errors, "mean_error");
  table.InsertColumn(3, v_estimates, "variance_estimate");
  table.InsertColumn(4, v_errors, "variance_error");
  //    table.InsertColumn(3, try_est, "error2");

  table.Save(output_path);
  return 0;
}