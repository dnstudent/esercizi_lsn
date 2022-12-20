#include <cmath>
#include <filesystem>
#include <iostream>
#include <iterator>
#include <optional>
#include <tuple>

#include <cxxopts.hpp>
#include <indicators/progress_bar.hpp>
#include <rapidcsv.h>

#include "ariel_random/ariel_random.hpp"
#include "config.hpp"
#include "molecular_systems/steppers/mc.hpp"
#include "molecular_systems/steppers/md.hpp"
#include "molecular_systems/system.hpp"
#include "options.hpp"
#include "utils.hpp"


#define SECTION "07"
#define EXERCISE SECTION "_4"

using std::string;
using Value = double;
static constexpr const bool tail_corrections = true;
static const std::array<std::string, 8> scalar_columns{"u_mean", "u_error", "e_mean", "e_error",
                                                       "T_mean", "T_error", "p_mean", "p_error"};
using MCSystem = LJMono<Value, tail_corrections, Ensamble::NVT>;
using MDSystem = LJMono<Value, tail_corrections, Ensamble::NVE>;
namespace ms_steppers = molecular_systems::steppers;
using namespace ex07;
using namespace indicators;

namespace co = cxxopts;
namespace fs = std::filesystem;


template<bool log_ar, class System, class Stepper>
void warmup(System &system, Method m, Ex4Options &p, Stepper &&stepper) {
    ms_steppers::StepSampler<Stepper> warmupper(std::forward<Stepper>(stepper));
    warmupper.sample(system, system.m_simulation.block_size * system.m_simulation.n_blocks);
    if constexpr (log_ar) { std::cout << warmupper.m_stepper.acceptance_rate(); }
    if (m == Method::MD) system.save_configurations(p.output_positions[m], p.output_velocities);
    else
        system.save_positions(p.output_positions[m]);
}

template<class System, class Stepper>
void take_measures(System &system, Method m, Ex4Options &p, Stepper &&stepper) {
    system.init_radial_func(p.n_bins);

    auto estimators = std::make_tuple(ProgAvg<Value>(), ProgAvg<Value>(), ProgAvg<Value>(),
                                      ProgAvg<Value>(), ProgAvg<std::vector<Value>>(p.n_bins));
    ms_steppers::BlockStats<Stepper, decltype(estimators), Variable::PotentialEnergy,
                            Variable::TotalEnergy, Variable::Temperature, Variable::Pressure,
                            Variable::RadialFn>
            block_stats(std::forward<Stepper>(stepper), std::move(estimators),
                        system.m_simulation.block_size);

    std::array<std::vector<Value>, scalar_columns.size()> scalar_results;
    for (size_t block = 0; block + 1 < system.m_simulation.n_blocks; block++) {
        const auto stats = utils::tuple_flatten(block_stats.statistics(system));
        utils::tuple_push_back(stats, scalar_results);
    }
    const auto stats = utils::tuple_flatten(block_stats.statistics(system));
    utils::tuple_push_back(stats, scalar_results);

    csv::Document scalar_table;
    utils::AppendColumns(scalar_table, scalar_columns, scalar_results);
    scalar_table.Save(p.output_dir[m] / "thermo.csv");

    csv::Document g_r_table;
    g_r_table.SetColumn(0, std::get<scalar_columns.size()>(stats));
    g_r_table.SetColumnName(0, "g_mean");
    g_r_table.SetColumn(1, std::get<scalar_columns.size() + 1>(stats));
    g_r_table.SetColumnName(1, "g_error");
    g_r_table.SetColumn(2, system.m_drs);
    g_r_table.SetColumnName(2, "r");
    g_r_table.Save(p.output_dir[m] / "g_r.csv");
}

int main(int argc, char const *argv[]) {
    // Defining some flags to control the program output
    cxxopts::Options options(EXERCISE, "How to run exercise " EXERCISE);
    // clang-format off
    options.add_options("Program")
      ("o,out", "Output dir", co::value<fs::path>()->default_value(RESULTS_DIR "/" SECTION "/"))
      ("in_mc", "Directory where data used to launch the run with the MC sampler is stored."
                "There must be at least a 'positions'"
                "file with the molecules positions and a 'settings' file with the simulation settings."
                "The presence of a 'velocities' file triggers a resume", co::value<fs::path>())
      ("in_md", "Directory where data used to launch the run with the MD integrator is stored."
                "There must be at least a 'positions'"
                "file with the molecules positions and a 'settings' file with the simulation settings."
                "The presence of a 'velocities' file triggers a resume", co::value<fs::path>())
      ("mc_settings", "Path to the MC settings file, if not present in 'in_mc'.", co::value<string>()->default_value(""))
      ("md_settings", "Path to the MD settings file, if not present in 'in_md'.", co::value<string>()->default_value(""))
      ("h,help", "Print this message");
    options.add_options("Rng seeding")
      ("p,primes_path", "Prime numbers path", co::value<string>()->default_value(PRIMES_PATH "Primes"))
      ("l,primes_line", "Line in primes_path to use", co::value<size_t>()->default_value("1"))
      ("s, seeds_path", "Seed path", co::value<string>()->default_value(""));
    options.add_options("Simulation")
      ("mc", "Whether to use the MC sampler", co::value<bool>())
      ("md", "Whether to use the MD sampler", co::value<bool>())
      ("warmup", "Whether it is a warmup or measure run", co::value<bool>())
      ("n,n_bins", "Number of bins for the radial function histogram", co::value<size_t>()->default_value("10"));
    // clang-format on
    auto user_params = options.parse(argc, argv);
    if (user_params.count("help")) {
        std::cout << options.help() << std::endl;
        exit(0);
    }

    // Define all the parameters
    const size_t PRIMES_LINE = user_params["l"].as<size_t>();
    const string PRIMES_SOURCE = user_params["p"].as<string>();
    Ex4Options p(user_params);
    const string SEEDS_SOURCE = p.rng_seed_path.string();

    auto rng = std::make_shared<ARandom>(SEEDS_SOURCE, PRIMES_SOURCE, PRIMES_LINE);
    auto m = Method::MC;
    if (p.sample[m]) {
        MCSystem mc_system(p.input_settings[m], p.input_positions[m]);
        auto stepper = ms_steppers::MC<Value, tail_corrections, ARandom>(
                mc_system.m_simulation.n_particles, mc_system.m_simulation.delta, rng);
        if (p.warmup) {
            warmup<true>(mc_system, m, p, std::move(stepper));
            //            std::cout << p.input_settings[m].string() << ": " << stepper.acceptance_rate()
            //                      << std::endl;
        } else {
            mc_system.init_velocities(*rng);
            take_measures(mc_system, m, p, std::move(stepper));
        }
        rng->SaveSeed((p.output_dir[m] / "rng.seed").string());
        mc_system.save_positions(p.output_positions[m]);
    }
    m = Method::MD;
    if (p.sample[m]) {
        MDSystem md_system(p.input_settings[m], p.input_positions[m]);
        if (p.resume[m]) {
            md_system.init_velocities(p.input_velocities);
        } else {
            md_system.init_velocities(*rng);
        }
        ms_steppers::MD2<Value, tail_corrections> stepper;
        if (p.warmup) {
            warmup<false>(md_system, m, p, std::move(stepper));
        } else {
            take_measures(md_system, m, p, std::move(stepper));
        }
        md_system.save_configurations(p.output_positions[m], p.output_velocities);
    }
    return 0;
}