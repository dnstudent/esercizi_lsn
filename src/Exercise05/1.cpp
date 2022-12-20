//
// Created by Davide Nicoli on 16/03/22.
//

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <iterator>
#include <numeric>
#include <optional>
#include <string>
#include <valarray>
#include <vector>

#include <cxxopts.hpp>
#include <rapidcsv.h>

#include "ariel_random/ariel_random.hpp"
#include "config.hpp"
#include "estimators/mean.hpp"
#include "samplers/MCMC/metropolis.hpp"
#include "structs.hpp"
#include "transitions/gauss.hpp"
#include "transitions/uniform.hpp"

using std::string;
using Value = double;
typedef double ProbSpace;
using StateSpace = std::valarray<Value>;
using samplers::mcmc::Metropolis;
using transitions::GaussNear;
using transitions::UniformNear;
namespace csv = rapidcsv;
namespace co = cxxopts;
namespace fs = std::filesystem;

#define SECTION "05"
#define EXERCISE SECTION "_1"

static const double LOG_COEFF_1S = -std::log(M_PI);
static const double LOG_COEFF_2P = -std::log(32 * M_PI);

inline Value radius(const StateSpace &x) {
    return std::sqrt(std::inner_product(std::begin(x), std::end(x), std::begin(x), double(0)));
}

// pdf: 1/π * e^-2r
struct psi_s_pdf {
    typedef ProbSpace prob_space;
    static prob_space logp(const StateSpace &x) { return -2 * radius(x) + LOG_COEFF_1S; }
};

// pdf: 1/32π * r^2 * e^-r * cos(th)^2
struct psi_2p_pdf {
    typedef ProbSpace prob_space;
    static prob_space logp(const StateSpace &x) {
        const auto r = radius(x);
        const auto costh = std::cos(std::atan2(x[1], x[0]));
        return LOG_COEFF_2P + 2 * std::log(r * std::abs(costh)) - r;
    }
};

/// Generate a sample of radii (atom distances from origin) using a specified probability distribution
/// \param sampler An MCMC sampler
/// \param rng A random number generator
/// \param sample_size The size of the sample
template<class MCMCSampler, class URBG>
void generate_estimates(MCMCSampler &sampler, URBG &rng, size_t sample_size, size_t n_blocks,
                        size_t warmup_steps, csv::Document &r_table) {
    std::vector<Value> radii(sample_size);
    std::vector<Value> r_estimates(n_blocks);
    std::vector<Value> r_errs(n_blocks);
    std::vector<double> acceptance_rates(n_blocks);
    estimators::ProgAvg<Value> estimator;
    sampler.warmup(warmup_steps, rng);

    for (size_t block = 0; block < n_blocks; block++) {
        acceptance_rates[block] = sampler.sample(radii.begin(), radii.end(), rng, radius);
        std::tie(r_estimates[block], r_errs[block]) = estimator(radii.cbegin(), radii.cend());
    }
    auto index = r_table.GetColumnCount();
    r_table.SetColumn(index, r_estimates);
    r_table.SetColumnName(index++, "mean");
    r_table.InsertColumn(index++, r_errs, "error");
    r_table.InsertColumn(index, acceptance_rates, "acceptance_rate");
}


template<class MCMCSampler, class URBG>
void generate_estimates_points(MCMCSampler &sampler, URBG &rng, size_t sample_size, size_t n_blocks,
                               size_t warmup_steps, csv::Document &r_table,
                               csv::Document &x_table) {
    std::vector<StateSpace> points(sample_size);
    std::vector<Value> radii(sample_size);
    std::vector<Value> r_estimates(n_blocks);
    std::vector<Value> r_errs(n_blocks);
    std::vector<double> acceptance_rates(n_blocks);
    estimators::ProgAvg<Value> estimator;
    sampler.warmup(warmup_steps, rng);

    for (size_t block = 0; block < n_blocks; block++) {
        acceptance_rates[block] = sampler.sample(points.begin(), points.end(), rng);
        std::transform(points.cbegin(), points.cend(), radii.begin(), radius);
        std::tie(r_estimates[block], r_errs[block]) = estimator(radii.cbegin(), radii.cend());
        for (auto &p: points) {
            x_table.InsertRow(x_table.GetRowCount(), std::vector(std::begin(p), std::end(p)));
        }
    }
    auto index = r_table.GetColumnCount();
    r_table.SetColumn(index, r_estimates);
    r_table.SetColumnName(index++, "mean");
    r_table.InsertColumn(index++, r_errs, "error");
    r_table.InsertColumn(index, acceptance_rates, "acceptance_rate");
}


template<class PDF, class StepSampler, class URBG>
void estimate_and_store(PDF &&pdf, const std::string_view pdf_name, StepSampler &&stepper,
                        const std::string_view stepper_name, const ExOptions<Value> &p, URBG &rng,
                        bool verbose) {
    if (verbose) std::cout << "Sampling " << pdf_name << ' ' << stepper_name << std::endl;
    Metropolis position_sampler(p.S0, std::forward<PDF>(pdf), std::forward<StepSampler>(stepper));
    csv::Document radii_table;
    const auto out_dir = p.output_dir / std::string(pdf_name) / std::string(stepper_name);
    if (!fs::exists(out_dir)) fs::create_directories(out_dir);
    if (p.save_positions) {
        csv::Document positions_table;
        positions_table.SetColumnName(0, "x");
        positions_table.SetColumnName(1, "y");
        positions_table.SetColumnName(2, "z");
        generate_estimates_points(position_sampler, rng, p.block_size, p.n_blocks, p.warmup_steps,
                                  radii_table, positions_table);
        positions_table.Save(out_dir / "positions.csv");
    } else {
        generate_estimates(position_sampler, rng, p.block_size, p.n_blocks, p.warmup_steps,
                           radii_table);
    }
    radii_table.Save(out_dir / "radii.csv");
}

int main(int argc, char const *argv[]) {
    // Defining some flags to control the program output
    cxxopts::Options options(EXERCISE, "How to run exercise " EXERCISE);
    // clang-format off
    options.add_options("Program")
      ("o,out", "Output directory", co::value<fs::path>()->default_value(RESULTS_DIR "/" SECTION "/"))
      ("h,help", "Print this message")
      ("v,verbose", "Whether to be verbose", co::value<bool>());
    options.add_options("Rng seeding")
      ("p,primes_path", "Prime numbers path", co::value<string>()->default_value(PRIMES_PATH "Primes"))
      ("l,primes_line", "Line in primes_path to use", co::value<size_t>()->default_value("1"))
      ("s, seeds_path", "Seed path", co::value<string>()->default_value(SEEDS_PATH "seed.in"));
    options.add_options("Sampling")
      ("M,n_throws", "Number of throws, i.e. total number of generated points", co::value<size_t>()->default_value("1000000"))
      ("N,n_blocks", "Number of blocks in which the throws are distributed", co::value<size_t>()->default_value("100"))
      ("w,n_warmup", "Number of MCMC warmup steps", co::value<size_t>()->default_value("1000"))
      ("steppers_config", "Four scalar values representing the steppers step sizes ordered as: s_uniform,s_gauss,2p_uniform,2p_gauss", co::value<std::vector<Value>>()->default_value("1.0,1.0,1.5,1.5"))
      ("starting_point", "Starting point for the MCMC sampler", co::value<std::vector<Value>>()->default_value("0.0,0.0,0.0"))
      ("uniform", "Perform the uniform sampling", co::value<bool>())
      ("gauss", "Perform the gaussian sampling", co::value<bool>())
      ("orbital_s", "Perform the sampling on orbital s", co::value<bool>())
      ("orbital_2p", "Perform the sampling on orbital 2p", co::value<bool>())
      ("positions", "Store sampled points in a csv file", co::value<bool>());
    // clang-format on
    auto user_params = options.parse(argc, argv);
    if (user_params.count("help")) {
        std::cout << options.help() << std::endl;
        exit(0);
    }
    const bool verbose = user_params["v"].as<bool>();

    // Defining all the parameters
    const string PRIMES_SOURCE = user_params["p"].as<string>();
    const size_t PRIMES_LINE = user_params["l"].as<size_t>();
    const string SEEDS_SOURCE = user_params["s"].as<string>();
    ExOptions<Value> p(user_params);


    // Instantiating the rng and MD integrator
    ARandom rng(SEEDS_SOURCE, PRIMES_SOURCE, PRIMES_LINE);

    const size_t d = p.S0.size();
    if (p.sample_s && p.sample_uniform)
        estimate_and_store(psi_s_pdf(), "orbital_s",
                           UniformNear<ProbSpace, StateSpace>(p.step_unif_s, d), "uniform", p, rng,
                           verbose);
    if (p.sample_s && p.sample_gauss)
        estimate_and_store(psi_s_pdf(), "orbital_s",
                           GaussNear<ProbSpace, StateSpace>(p.step_gauss_s, d), "gauss", p, rng,
                           verbose);
    if (p.sample_2p && p.sample_uniform)
        estimate_and_store(psi_2p_pdf(), "orbital_2p",
                           UniformNear<ProbSpace, StateSpace>(p.step_unif_2p, d), "uniform", p, rng,
                           verbose);
    if (p.sample_2p && p.sample_gauss)
        estimate_and_store(psi_2p_pdf(), "orbital_2p",
                           GaussNear<ProbSpace, StateSpace>(p.step_gauss_2p, d), "gauss", p, rng,
                           verbose);

    return 0;
}