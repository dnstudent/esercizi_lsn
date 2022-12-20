#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <iterator>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include <cxxopts.hpp>
#include <rapidcsv.h>

#include "ariel_random/ariel_random.hpp"
#include "config.hpp"
#include "models/ising/1D/ising.hpp"
#include "models/ising/1D/simulation.hpp"
#include "models/ising/1D/variables.hpp"
#include "models/ising/structs.hpp"
#include "samplers/MCMC/gibbs.hpp"
#include "samplers/MCMC/metropolis.hpp"
#include "structs.hpp"
#include "transitions/uniform.hpp"

using std::string;
namespace csv = rapidcsv;
namespace co = cxxopts;
namespace fs = std::filesystem;
using namespace ising;
using namespace samplers::mcmc;
typedef double ProbSpace;
typedef double VarSpace;

// Variable types, using aliases
using Ising1D = Ising<VarSpace, 1>;

#define SECTION "06"
#define EXERCISE SECTION "_equilibrate"

/**
 * Utility function to launch the equilibration. Two runs are performed to cope with the possibility
 * of running in a local minimum and have biased results.
 * Check out the D1::Equilibrator documentation for additional insights.
 * @tparam compute_H Flag controlling wether to compute H.
 * @tparam compute_s Flag controlling wether to compute Sum_s.
 * @tparam compute_s2 Flag controlling wether to compute Sum_s2.
 * @tparam Sampler Sampler class.
 * @tparam URBG Random number generator class.
 * @param p Program options.
 * @param h External magnetic field.
 * @param sampler_name Tag used to identify runs.
 * @param rng1 Random number generator for the first sequence.
 * @param rng2 Random number generator for the second sequence.
 */
template<bool compute_H, bool compute_s, bool compute_s2, class Sampler, class URBG>
void equilibrate(ExOptions<VarSpace> &p, VarSpace h, std::string_view sampler_name, URBG &rng1,
                 URBG &rng2) {
    // First run
    auto ising_1 = std::make_shared<Ising1D>(p.n_spins, rng1, p.J, h, p.T);
    D1::Equilibrator<compute_H, compute_s, compute_s2, VarSpace, Sampler> equilibrator1(p.n_steps,
                                                                                        ising_1);
    equilibrator1.run(p.warmup_steps, rng1);
    equilibrator1.save_results(
            p.output_dir / (std::string(sampler_name) + "_" + std::to_string(h) + "_warmup1.csv"));
    if (p.save_spins)
        equilibrator1.save_state(p.output_dir / (std::string(sampler_name) + "_" +
                                                 std::to_string(h) + "_spins.csv"));
    // Second run
    auto ising_2 = std::make_shared<Ising1D>(p.n_spins, rng2, p.J, h, p.T);
    D1::Equilibrator<compute_H, compute_s, compute_s2, VarSpace, Sampler> equilibrator2(p.n_steps,
                                                                                        ising_2);
    equilibrator2.run(p.warmup_steps, rng2);
    equilibrator2.save_results(
            p.output_dir / (std::string(sampler_name) + "_" + std::to_string(h) + "_warmup2.csv"));
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
    options.add_options("Ising")
      ("n_spins", "Number of spins", co::value<size_t>()->default_value("50"))
      ("J,coupling", "Spin coupling in reduced units", co::value<VarSpace>()->default_value("1"))
      ("B,external_field", "External magnetic field in reduced units", co::value<VarSpace>()->default_value("0"))
      ("T,temperature", "Temperature in reduced units", co::value<VarSpace>()->default_value("0.5"));
    options.add_options("Sampling")
      ("M,n_steps", "Number of throws, i.e. total number of generated points", co::value<size_t>())
      ("S,block_size", "Size of the blocks in which the throws are distributed", co::value<size_t>()->default_value("1"))
      ("w,n_warmup", "Number of MCMC warmup steps", co::value<size_t>())
      ("metropolis", "Whether to sample using the Metropolis algorithm", co::value<bool>())
      ("gibbs", "Whether to sampler using the Gibbs algorithm", co::value<bool>())
      ("save_spins", "Whether to save the state of the Ising model", co::value<bool>())
      ("resume", "Whether to resume a previous run in the chosen output directory", co::value<bool>());
    // clang-format on
    auto user_params = options.parse(argc, argv);
    if (user_params.count("help")) {
        std::cout << options.help() << std::endl;
        exit(0);
    }

    // Defining all the parameters
    const string PRIMES_SOURCE = user_params["p"].as<string>();
    const size_t PRIMES_LINE = user_params["l"].as<size_t>();
    const string SEEDS_SOURCE = user_params["s"].as<string>();
    ExOptions<VarSpace> p(user_params);

    ARandom rng1(SEEDS_SOURCE, PRIMES_SOURCE, PRIMES_LINE);
    ARandom rng2(SEEDS_SOURCE, PRIMES_SOURCE, PRIMES_LINE + 1);


    if (p.metropolis) {
        equilibrate<true, false, true, SystemMetropolis<VarSpace>>(p, 0., "metropolis", rng1, rng2);
        equilibrate<false, true, false, SystemMetropolis<VarSpace>>(p, p.h, "metropolis", rng1,
                                                                    rng2);
    }
    if (p.gibbs) {
        equilibrate<true, false, true, SystemGibbs<VarSpace>>(p, 0., "gibbs", rng1, rng2);
        equilibrate<false, true, false, SystemGibbs<VarSpace>>(p, p.h, "gibbs", rng1, rng2);
    }
    return 0;
}