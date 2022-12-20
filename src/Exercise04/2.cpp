#include <cmath>
#include <filesystem>
#include <iostream>
#include <iterator>
#include <optional>

#include <cxxopts.hpp>
#include <indicators/progress_bar.hpp>
#include <rapidcsv.h>

#include "ariel_random/ariel_random.hpp"
#include "config.hpp"
#include "molecular_systems/steppers/md.hpp"
#include "molecular_systems/system.hpp"


#define SECTION "04"
#define EXERCISE SECTION "_2"

#define N_VARS 5

using std::string;
using namespace indicators;
using Value = double;
using namespace molecular_systems::steppers;
namespace co = cxxopts;
namespace fs = std::filesystem;

int main(int argc, char const *argv[]) {
    // Defining some flags to control the program output
    cxxopts::Options options(EXERCISE, "How to run exercise " EXERCISE);
    // clang-format off
    options.add_options("Program")
      ("o,out", "Output dir", co::value<fs::path>()->default_value(RESULTS_DIR "/" SECTION "/"))
      ("h,help", "Print this message");
    options.add_options("Rng seeding")
      ("p,primes_path", "Prime numbers path", co::value<string>()->default_value(PRIMES_PATH "Primes"))
      ("l,primes_line", "Line in primes_path to use", co::value<size_t>()->default_value("1"))
      ("s, seeds_path", "Seed path", co::value<string>()->default_value(SEEDS_PATH "seed.in"));
    options.add_options("Simulation")
      ("settings", "Path to the simulator settings", co::value<fs::path>()->default_value(MD_SETTINGS_PATH "input.solid"))
      ("x,configuration", "Path to molecular configuration. Must be three columns representing each molecule's position components", co::value<fs::path>()->default_value(LATTICES_PATH "config.fcc"))
      ("v,velocities", "Path to molecular velocities, to resume transform previous run", co::value<std::string>()->default_value(""))
      ("N,save_every", "Save every N frames", co::value<size_t>()->default_value("0"));
    // clang-format on
    auto user_params = options.parse(argc, argv);
    if (user_params.count("help")) {
        std::cout << options.help() << std::endl;
        exit(0);
    }

    // Defines all the parameters
    const auto CONFIGURATION_PATH = user_params["x"].as<fs::path>();
    const auto SETTINGS_PATH = user_params["settings"].as<fs::path>();
    const auto VELOCITIES_STRING = user_params["v"].as<std::string>();
    std::optional<fs::path> VELOCITIES_PATH;
    if (VELOCITIES_STRING.empty()) {
        VELOCITIES_PATH = std::nullopt;
    } else {
        VELOCITIES_PATH = fs::path(VELOCITIES_STRING);
    }
    const auto SAVE_EVERY_N_FRAMES = user_params["N"].as<size_t>();
    const auto OUTPUT_DIR = user_params["o"].as<fs::path>();
    if (!fs::exists(OUTPUT_DIR)) { fs::create_directories(OUTPUT_DIR); }
    const string PRIMES_SOURCE = user_params["p"].as<string>();
    const size_t PRIMES_LINE = user_params["l"].as<size_t>();
    const string SEEDS_SOURCE = user_params["s"].as<string>();
    const auto FRAMES_DIR = OUTPUT_DIR / "frames/";
    if (!fs::exists(FRAMES_DIR)) fs::create_directories(FRAMES_DIR);

    // Instantiating the rng and MD integrator
    ARandom rng(SEEDS_SOURCE, PRIMES_SOURCE, PRIMES_LINE);
    auto system = std::make_shared<LJMono<Value, false, Ensamble::NVE>>(
            SETTINGS_PATH, CONFIGURATION_PATH, VELOCITIES_PATH, rng);
    MD integrator(system);

    // Vectors where mean estimations and their variances will be stored
    std::array<std::vector<Value>, 3 * N_VARS> stats;
    // A progressbar
    ProgressBar pbar(option::MaxProgress{system->m_simulation.n_blocks},
                     option::ShowElapsedTime{true}, option::ShowRemainingTime{true},
                     option::BarWidth{80});
    // Computing and storing block statistics
    for (size_t i = 0UL; i < system->m_simulation.n_blocks; i++) {
        const auto block_results = integrator.block_estimates(SAVE_EVERY_N_FRAMES, FRAMES_DIR);
        for (size_t var = 0UL; var < N_VARS; var++) {
            stats[3 * var].push_back(std::get<0>(block_results[var]));
            stats[3 * var + 1].push_back(std::get<1>(block_results[var]));
            stats[3 * var + 2].push_back(std::get<2>(block_results[var]));
        }
        pbar.tick();
    }
    // Saving final molecular positions and velocities
    system->save_configurations((OUTPUT_DIR / "config.positions"),
                                (OUTPUT_DIR / "config.velocities"));

    // Storing results in a csv file
    const auto variable_names = LJMono<Value, false, Ensamble::NVE>::variable_names();
    rapidcsv::Document table;
    for (size_t var = 0UL; var < N_VARS; var++) {
        table.InsertColumn(3 * var, stats[3 * var], variable_names[var] + "_blockmean");
        table.InsertColumn(3 * var + 1, stats[3 * var + 1], variable_names[var] + "_progmean");
        table.InsertColumn(3 * var + 2, stats[3 * var + 2], variable_names[var] + "_error");
    }
    table.RemoveColumn(table.GetColumnCount() - 1);
    table.Save(OUTPUT_DIR / "thermo.csv");
    rng.SaveSeed((OUTPUT_DIR / "rng.seed").string());
    return 0;
}