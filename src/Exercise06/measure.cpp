//
// Created by Davide Nicoli on 16/03/22.
//

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iterator>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include <cxxopts.hpp>

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
using U_var = D1::InternalEnergy<VarSpace, ProgAvg<VarSpace>>;
using C_var = D1::HeatCapacity<VarSpace, ProgVariance<VarSpace>>;
using X_var = D1::MagneticSusceptivity<VarSpace, ProgAvg<VarSpace>>;
using M_var = D1::Magnetization<VarSpace, ProgAvg<VarSpace>>;
using Ising1D = Ising<VarSpace, 1>;

#define SECTION "06"
#define EXERCISE SECTION "_1"

/**
 * Path to the last sampled state
 * @param p Program options
 * @param h Value of the external magnetic field.
 * @param sampler_name The name of the sampler which was used.
 * @return The path to the last sampled state.
 */
fs::path state_path(ExOptions<VarSpace> &p, VarSpace h, std::string_view sampler_name) {
    return p.output_dir / (std::string(sampler_name) + "_" + std::to_string(h) + "_spins.csv");
}
/**
 * A helper function to perform measures. The measures results will be stored in csv files named
 * "{sampler_name}_{h}_variables.csv"
 * @tparam Sampler The class of the MCMC algorithm which will be used.
 * @tparam URBG The class of the random number generator.
 * @tparam compute_H Flag to enable/disable the computation of H at compile time.
 * @tparam compute_s Flag to enable/disable the computation of the sum of spins at compile time.
 * @tparam compute_s2 Flag to enable/disable the computation of the sum of the product of spins at compile time.
 * @tparam ThermoVars Classes of the thermodynamic variables.
 * @param p Program options.
 * @param h Value of the external magnetic field.
 * @param sampler_name A tag used to identify the run.
 * @param rng The random number generator.
 */
template<class Sampler, class URBG, bool compute_H, bool compute_s, bool compute_s2,
         class... ThermoVars>
void measure(ExOptions<VarSpace> &p, VarSpace h, std::string_view sampler_name, URBG &rng) {
    // Initialization of the Ising model: are we resuming a previous run?
    std::shared_ptr<Ising1D> ising_model;
    if (p.resume)
        ising_model = std::make_shared<Ising1D>(state_path(p, h, sampler_name), p.J, h, p.T);
    else
        ising_model = std::make_shared<Ising1D>(p.n_spins, rng, p.J, h, p.T);
    // Initialization of the simulator. It is instructed on which variables must be measured (ThermoVars) and
    // which intermediate variables should be computed (compute_H,...) for performance reasons
    D1::Simulator<compute_H, compute_s, compute_s2, VarSpace, Sampler, ThermoVars...> simulator(
            p.block_size, ising_model, ThermoVars(*ising_model)...);
    simulator.run(p.n_blocks, p.warmup_steps, rng);
    simulator.save_results(p.output_dir / (std::string(sampler_name) + "_" + std::to_string(h) +
                                           "_variables.csv"));
    if (p.save_spins) ising_model->save_state(state_path(p, h, sampler_name));
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
      ("B,external_field", "External magnetic field in reduced units (a value != 0 will be used only to compute M)", co::value<VarSpace>()->default_value("0.02"))
      ("T,temperature", "Temperature in reduced units", co::value<VarSpace>()->default_value("0.5"))
      ("resume_from", "State from which sampling will be started", co::value<fs::path>());
    options.add_options("Sampling")
      ("M,n_steps", "Number of MC steps, i.e. number of blocks", co::value<size_t>())
      ("S,block_size", "Size of the blocks in which the throws are distributed", co::value<size_t>())
      ("w,n_warmup", "Number of MCMC warmup steps", co::value<size_t>()->default_value("0"))
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

    // Storing the flag values
    const string PRIMES_SOURCE = user_params["p"].as<string>();
    const size_t PRIMES_LINE = user_params["l"].as<size_t>();
    const string SEEDS_SOURCE = user_params["s"].as<string>();
    ExOptions<VarSpace> p(user_params);

    ARandom rng(SEEDS_SOURCE, PRIMES_SOURCE, PRIMES_LINE);
    if (p.metropolis) {
        // The case in which the external field is globally set to 0 must be distinguished, otherwise
        // the output files will overwrite each other.
        if (p.h == 0)
            measure<SystemMetropolis<VarSpace>, decltype(rng), true, true, true, U_var, C_var,
                    X_var, M_var>(p, 0., "metropolis", rng);
        else {
            // The computation of U, C and X only requires the evaluation of H and Sum_s2, hence Sum_s is disabled
            measure<SystemMetropolis<VarSpace>, decltype(rng), true, false, true, U_var, C_var,
                    X_var>(p, 0., "metropolis", rng);
            // The computation of M only requires the evaluation of Sum_s, hence H and Sum_s2 are disabled
            measure<SystemMetropolis<VarSpace>, decltype(rng), false, true, false, M_var>(
                    p, p.h, "metropolis", rng);
        }
    }
    if (p.gibbs) {
        if (p.h == 0)
            measure<SystemGibbs<VarSpace>, decltype(rng), true, true, true, U_var, C_var, X_var,
                    M_var>(p, 0., "gibbs", rng);
        else {
            measure<SystemGibbs<VarSpace>, decltype(rng), true, false, true, U_var, C_var, X_var>(
                    p, 0., "gibbs", rng);
            measure<SystemGibbs<VarSpace>, decltype(rng), false, true, false, M_var>(p, p.h,
                                                                                     "gibbs", rng);
        }
    }
    return 0;
}