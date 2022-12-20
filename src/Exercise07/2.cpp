#include <cmath>
#include <filesystem>
#include <iostream>
#include <iterator>
#include <optional>

#include <cxxopts.hpp>
#include <rapidcsv.h>

#include "ariel_random/ariel_random.hpp"
#include "config.hpp"
#include "molecular_systems/steppers/collectors.hpp"
#include "molecular_systems/steppers/mc.hpp"
#include "molecular_systems/system.hpp"
#include "options.hpp"
#include "utils.hpp"

#define SECTION "07"
#define EXERCISE SECTION "_2"

#define VARIABLES                                                                                  \
    Variable::PotentialEnergy, Variable::KineticEnergy, Variable::TotalEnergy,                     \
            Variable::Temperature, Variable::Pressure
using std::string;
using Value = double;
static constexpr bool tail_corrections = true;
using MCSystem = LJMono<Value, tail_corrections, Ensamble::NVT>;
namespace ms_step = molecular_systems::steppers;
using namespace ex07;
using Outs = MeasureOutputs<Value, VARIABLES>;

namespace co = cxxopts;
namespace fs = std::filesystem;

int main(int argc, char const *argv[]) {
    // Defining some flags to control the program i/o
    cxxopts::Options options(EXERCISE, "How to run exercise " EXERCISE);
    // clang-format off
    options.add_options("Program")
      ("o,out", "Output dir", co::value<fs::path>()->default_value(RESULTS_DIR "/" SECTION "/"))
      ("i,in", "Directory where data used to launch the run is stored. There must be at least a 'positions'"
               "file with the molecules positions.", co::value<fs::path>())
      ("settings", "Path to the settings file, if not present in 'in'.", co::value<string>()->default_value(""))
      ("h,help", "Print this message");
    options.add_options("Rng seeding")
      ("p,primes_path", "Prime numbers path", co::value<string>()->default_value(PRIMES_PATH "Primes"))
      ("l,primes_line", "Line in primes_path to use", co::value<size_t>()->default_value("0"))
      ("s,seeds_path", "Seed path", co::value<string>()->default_value(""));
    options.add_options("Simulation")
      ("N,save_every", "Save every N frames", co::value<size_t>()->default_value("0"));
    // clang-format on
    auto user_params = options.parse(argc, argv);
    if (user_params.count("help")) {
        std::cout << options.help() << std::endl;
        exit(0);
    }

    const size_t PRIMES_LINE = user_params["l"].as<size_t>();
    const string PRIMES_SOURCE = user_params["p"].as<string>();
    Ex2Options p(user_params);

    // Initializing the molecular system
    auto rng = std::make_shared<ARandom>(p.rng_seed_path.string(), PRIMES_SOURCE, PRIMES_LINE);

    MCSystem system(p.settings_path, p.positions_path);
    // Initializing the stepper, which evolves the system given a certain logic
    ms_step::MC<Value, tail_corrections, ARandom> stepper(system.m_simulation.n_particles,
                                                          system.m_simulation.delta, rng);
    // Initializing the sampler, which performes instantaneous measurements while evolving the system
    ms_step::StepSampler<decltype(stepper), VARIABLES> sampler(std::move(stepper));
    // Performing the measurements
    const Outs measures = sampler.sample(system, system.m_simulation.block_size);

    // Saving the final configuration...
    system.save_positions(p.output_positions);
    // ... random seed ...
    rng->SaveSeed((p.output_dir / "rng.seed").string());
    std::cout << sampler.m_stepper.acceptance_rate();
    // ... and measurements results
    csv::Document table;
    table.InsertColumn(0, measures.get_measures<PotentialEnergy>(), "U/N");
    table.InsertColumn(1, measures.get_measures<KineticEnergy>(), "K/N");
    table.InsertColumn(2, measures.get_measures<TotalEnergy>(), "E/N");
    table.InsertColumn(3, measures.get_measures<Temperature>(), "T");
    table.InsertColumn(4, measures.get_measures<Pressure>(), "p");
    table.RemoveColumn(table.GetColumnCount() - 1);
    table.Save(p.output_dir / "thermo.csv");
    return 0;
}