//
// Created by Davide Nicoli on 25/10/22.
//

#include <memory>
#include <utility>
#include <valarray>
#include <vector>

#include <cxxopts.hpp>
#include <indicators/progress_bar.hpp>
#include <rapidcsv.h>

#include "config.hpp"
#include "distributions/exercises.hpp"
#include "estimators/estimators.hpp"
#include "mc_integrators/integrator.hpp"
#include "options.hpp"
#include "transitions/gauss.hpp"
#include "transitions/uniform.hpp"
#include "utils.hpp"
#include "variational_mc/simulated_annealing.hpp"

#define SECTION "08"
#define EXERCISE SECTION "_2"

using std::string;
using field = double;
using param_space = std::array<field, 2>;
using prob_space = double;
namespace co = cxxopts;
namespace fs = std::filesystem;

using namespace ex08;
using namespace variational_mc;
using namespace mc_integrator;
using namespace distributions::ex08;
using namespace functions::ex08;
using namespace transitions;


void set_options(cxxopts::Options &options) {
    // clang-format off
    options.add_options("Program")
      ("o,out", "Path to the directory in which outputs will be stored", co::value<fs::path>()->default_value(RESULTS_DIR "/" SECTION))
      ("h,help", "Print this message");
    options.add_options("Annealing")
      ("N,n_blocks", "Number of blocks for the estimation of <H>", co::value<size_t>()->default_value("100"))
      ("W,block_size", "Size of the blocks", co::value<size_t>()->default_value("10000"))
      ("b,bounds", "Bounds for the histogram of ψ", co::value<std::vector<field>>()->default_value("-3,3"))
      ("B,n_bins", "Number of bins used to produce the histogram of ψ", co::value<size_t>()->default_value("100"))
      ("mu", "ψ's µ", co::value<field>())
      ("sigma", "ψ's σ", co::value<field>());
    options.add_options("Rng seeding")
      ("p,primes_path", "Prime numbers path", co::value<string>()->default_value(PRIMES_PATH "Primes"))
      ("l,primes_line", "Line in primes_path to use", co::value<size_t>()->default_value("1"))
      ("s, seeds_path", "Seed path", co::value<string>()->default_value(SEEDS_PATH "seed.in"));
    // clang-format on
}

template<class URBG>
auto compute_estimate(const ExPsiOptions<field> &p, URBG &rng) {
    Integrand<field> Hpsi(p.mu, p.sigma);
    Metropolis psi_sampler(/*start=*/field(0), /*pdf=*/Trial<field>(p.mu, p.sigma),
                           /*transition=*/Uniform<double, field>(p.a, p.b));
    std::vector<field> psi_x_sample(p.n_blocks * p.block_size);
    std::vector<field> estimate(p.n_blocks);
    std::vector<field> error(estimate.size());
    {
        Integrator<field, decltype(psi_sampler)> I(std::move(psi_sampler));
        I(Hpsi, estimate.begin(), estimate.end(), error.begin(), psi_x_sample.begin(),
          psi_x_sample.end(), rng);
    }
    {
        csv::Document table;
        utils::AppendColumns(table, {"H_estimate", "H_error"}, std::make_tuple(estimate, error));
        table.Save(p.out / "H_min.csv");
    }
    {
        std::vector<size_t> bins(p.n_bins);
        std::vector<field> bins_edges(p.n_bins);
        utils::histogram(psi_x_sample.cbegin(), psi_x_sample.cend(), bins.begin(), bins.end(),
                         bins_edges.begin(), p.a, p.b);
        csv::Document table;
        utils::AppendColumns(table, {"psi", "l_edge"}, std::make_tuple(bins, bins_edges));
        table.Save(p.out / "psi.csv");
    }
}


int main(int argc, char const *argv[]) {
    // Defining some flags to control the program output
    cxxopts::Options options(EXERCISE, "How to run exercise " EXERCISE);
    set_options(options);
    auto user_params = options.parse(argc, argv);
    if (user_params.count("help")) {
        std::cout << options.help() << std::endl;
        exit(0);
    }

    // Define all the parameters
    const size_t PRIMES_LINE = user_params["l"].as<size_t>();
    const string PRIMES_SOURCE = user_params["p"].as<string>();
    const string SEEDS_SOURCE = user_params["s"].as<string>();
    ExPsiOptions<field> p(user_params);
    ARandom rng(SEEDS_SOURCE, PRIMES_SOURCE, PRIMES_LINE);
    compute_estimate(p, rng);

    return 0;
}