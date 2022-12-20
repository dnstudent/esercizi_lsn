//
// Created by Davide Nicoli on 25/10/22.
//
#include <array>
#include <iostream>
#include <memory>
#include <utility>
#include <vector>

#include <cxxopts.hpp>
#include <indicators/progress_bar.hpp>
#include <rapidcsv.h>

#include "config.hpp"
#include "distributions/exercises.hpp"
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

template<class URBG>
class Loss {
public:
    Loss(size_t n_blocks, size_t block_size, std::shared_ptr<URBG> rng)
        : m_n_blocks(n_blocks), m_rng(std::move(rng)), m_block(block_size), m_transf(block_size) {}

    auto operator()(const param_space &p) {
        if (p[0] < 0 || p[1] <= 0)
            return std::make_tuple(std::numeric_limits<field>::infinity(),
                                   std::numeric_limits<field>::quiet_NaN());
        m_H_estimator.reset();
        Metropolis sampler(/*start=*/static_cast<field>(p[0]),
                           /*pdf=*/Trial<field>(/*mu=*/p[0], /*sigma=*/p[1]),
                           // this value provides an acceptance rate of ~0.5
                           /*transition=*/UniformNear<prob_space, field>(2.5));
        // Warmup is not needed as the sampler already starts from the maximum probability state
        Integrand<field> Hpsi(/*mu=*/p[0], /*sigma=*/p[1]);
        for (size_t i = 0; i + 1 < m_n_blocks; i++) {
            sampler.sample(m_block.begin(), m_block.end(), *m_rng);
            std::transform(m_block.cbegin(), m_block.cend(), m_transf.begin(), Hpsi);
            m_H_estimator(m_transf.cbegin(), m_transf.cend());
        }
        sampler.sample(m_block.begin(), m_block.end(), *m_rng);
        std::transform(m_block.cbegin(), m_block.cend(), m_transf.begin(), Hpsi);
        return m_H_estimator(m_transf.cbegin(), m_transf.cend());
    }

private:
    const size_t m_n_blocks;
    std::shared_ptr<URBG> m_rng;
    std::vector<field> m_block;
    std::vector<field> m_transf;
    ProgAvg<field> m_H_estimator{};
};

void set_options(cxxopts::Options &options) {
    // clang-format off
    options.add_options("Program")
      ("o,out", "Path to the directory in which outputs will be stored", co::value<fs::path>()->default_value(RESULTS_DIR "/" SECTION))
      ("h,help", "Print this message");
    options.add_options("Annealing")
      ("N,n_steps", "Number of temperature variations", co::value<size_t>()->default_value("10"))
      ("W,n_explore", "Number of steps taken while keeping the temperature constant", co::value<size_t>()->default_value("10"))
      ("n_blocks", "Number of blocks used to estimate <H> in a single exploration step", co::value<size_t>()->default_value("100"))
      ("block_size", "Block size used to estimate <H> in a single exploration step", co::value<size_t>()->default_value("100"))
      ("T0", "Starting temperature", co::value<field>()->default_value("10"))
      ("Tf", "Final temperature", co::value<field>()->default_value("0.0001"))
      ("p0", "Initial parameter guess", co::value<std::vector<field>>()->default_value("1,1"))
      ("stddev", "Stddev of the normal distribution used to sample the next pair of parameters", co::value<field>()->default_value("0.05"));
    options.add_options("Rng seeding")
      ("p,primes_path", "Prime numbers path", co::value<string>()->default_value(PRIMES_PATH "Primes"))
      ("l,primes_line", "Line in primes_path to use", co::value<size_t>()->default_value("1"))
      ("s, seeds_path", "Seed path", co::value<string>()->default_value(SEEDS_PATH "seed.in"));
    // clang-format on
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
    ExOptions<field> p(user_params);

    auto rng = std::make_shared<ARandom>(SEEDS_SOURCE, PRIMES_SOURCE, PRIMES_LINE);
    ProgAvg<field> H_estimator;
    SimulatedAnnealing sa(
            /*loss_fn=*/Loss(p.n_blocks, p.block_size, rng),
            /*q=*/GaussNear<prob_space, param_space>(p.stddev));
    const size_t trajectory_size = p.n_T_steps * p.n_explore_steps + 1;
    std::vector<param_space> params(trajectory_size);
    std::vector<field> energies(trajectory_size);
    std::vector<field> energies_errors(trajectory_size);
    std::vector<field> temperatures(trajectory_size);
    LogScheduler<field> T_scheduler(p.T0, p.Tf, p.n_T_steps);
    sa.anneal(/*p0=*/{p.m0, p.s0}, p.n_explore_steps, T_scheduler, params.begin(), energies.begin(),
              energies_errors.begin(), temperatures.begin(), *rng);


    std::vector<field> mus(trajectory_size);
    std::vector<field> sigmas(trajectory_size);
    std::transform(params.cbegin(), params.cend(), mus.begin(),
                   [](const auto param) { return param[0]; });
    std::transform(params.cbegin(), params.cend(), sigmas.begin(),
                   [](const auto param) { return param[1]; });

    csv::Document params_table;
    utils::AppendColumns(params_table, {"H_estimate", "H_error", "mu", "sigma", "T"},
                         std::make_tuple(energies, energies_errors, mus, sigmas, temperatures));
    params_table.Save(p.out / "annealing.csv");
    return 0;
}