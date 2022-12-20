#include <cmath>
#include <exception>
#include <fstream>
#include <memory>
#include <numeric>
#include <random>
#include <valarray>
#include <vector>

#include <cxxopts.hpp>
#include <indicators/progress_bar.hpp>

#include "ariel_random/ariel_random.hpp"
#include "config.hpp"
#include "genetic_process.hpp"
#include "genetic_tsp/crossovers.hpp"
#include "genetic_tsp/tsp_ga.hpp"
#include "options.hpp"
#include "utils.hpp"

namespace co = cxxopts;
using std::string;
using namespace indicators;

#define SECTION "09"
#define EXERCISE SECTION "_1"
#define N_CITIES 34

namespace csv = rapidcsv;


void set_options(cxxopts::Options &options) {
    // clang-format off
    options.add_options("Program")
      ("o,out", "Path to the directory in which outputs will be stored", co::value<fs::path>()->default_value(RESULTS_DIR "/" SECTION))
      ("i,in", "Path to the csv file of coordinates", co::value<fs::path>())
      ("crossover", "Which crossover algorithm(s) to use: 'ex' for the one proposed with the exercises, 'exmod' for that but modified, "
       "'my1' or 'my2' for my algorithms.", co::value<string>())
      ("h,help", "Print this message");
    options.add_options("Genetic")
      ("n,n_iter", "Number of iterations", co::value<size_t>()->default_value("5000"))
      ("m,pop_size", "Population size", co::value<size_t>()->default_value("40000"))
      ("r,mut_rate", "Mutation rate", co::value<double>()->default_value("0.05"))
      ("f,fusion_p", "Fusion algorithm's probability of choosing my2 vs exmod", co::value<double>()->default_value("0.7"));
    options.add_options("Rng seeding")
      ("p,primes_path", "Prime numbers path", co::value<string>()->default_value(PRIMES_PATH "Primes"))
      ("l,primes_line", "Line in primes_path to use", co::value<size_t>()->default_value("0"))
      ("s, seeds_path", "Seed path", co::value<string>()->default_value(SEEDS_PATH "seed.in"));
    // clang-format on
}


template<typename Coordinates, class Crossover, class URBG>
auto make_run_gp(Coordinates &&coordinates, Crossover &&crossover, ex09::ExOptions &p, URBG &rng) {
    // Initializing the genetic process with the algorithm "TSP" (Travelling Salesman Problem)
    // The Process object implements the loops that evolve the system, while TSP describes the
    // individuals and the operations they can undergo.
    genetic::Process gp(
            (TSP(std::forward<Coordinates>(coordinates), std::forward<Crossover>(crossover))));

    using Individual = typename decltype(gp)::Individual;

    // Defining two containers for the population and its fitnesses.
    std::vector<Individual> population(p.pop_size);
    std::vector<double> evaluations(p.pop_size);
    std::vector<double> distances(p.n_iter);
    ProgressBar pbar{option::MaxProgress{p.n_iter}, option::ShowElapsedTime{true},
                     option::ShowRemainingTime{true}, option::BarWidth{80},
                     option::PrefixText{"cross_tsp: " + ex09::tag_from(p.algo)}};
    // Running the genetic process; the results are stored in population and evaluations
    gp.template run<true>(population.begin(), p.pop_size, evaluations.begin(), p.n_iter, p.mut_rate,
                          distances.begin(), pbar, rng);
    auto [best_found, best_fitness] = gp.get_best();
    population.push_back(best_found);
    evaluations.push_back(best_fitness);
    return std::make_tuple(population, evaluations, coordinates, distances);
}

auto generate_and_run_gp(ex09::CrossAlgo algo, ex09::ExOptions &p) {
    using point = std::array<double, 2>;

    ARandom rng(p.seeds_path, p.primes_path, p.primes_line);
    // Generating the coordinates for the circle problem
    std::array<point, N_CITIES> coordinates{};
    genetic::load_coordinates(p.in_path, coordinates, /*header=*/false);
#define RUN_GP(crossover) make_run_gp(std::move(coordinates), crossover, p, rng)
    if (algo == ex09::Exercise) return RUN_GP(ExerciseCrossover<N_CITIES>{});
    else if (algo == ex09::ExerciseMod)
        return RUN_GP(ExerciseModCrossover<N_CITIES>{});
    else if (algo == ex09::MyAlgo2)
        return RUN_GP((MyCrossover2<N_CITIES>{}));
    else if (algo == ex09::Fusion)
        return RUN_GP((Fusion<N_CITIES>(p.fusion_p)));
    else if (algo == ex09::Dummy)
        return RUN_GP(CloneCrossover{});
    throw std::runtime_error("Not implemented");
}

int main(int argc, char *argv[]) {
    // Defining some flags to control the program output
    cxxopts::Options options(EXERCISE, "How to run exercise " EXERCISE);
    set_options(options);
    auto user_params = options.parse(argc, argv);
    if (user_params.count("help")) {
        std::cout << options.help() << std::endl;
        exit(0);
    }

    // A struct encoding the options given to the program
    ex09::ExOptions p(user_params);


    auto [population, evaluations, coordinates, distances] = generate_and_run_gp(p.algo, p);
    csv::Document table;
    for (size_t i = 0; i < population.size(); i++) {
        auto row = std::vector<size_t>(population[i].begin(), population[i].end());
        row.insert(row.cbegin(), 0UL);
        table.SetRow(i, row);
    }
    std::for_each(evaluations.begin(), evaluations.end(), [](auto &val) { val = 1.0 / val; });
    utils::AppendColumns(table, {"total_distance"}, std::make_tuple(evaluations));
    table.Save(p.out_dir / (ex09::tag_from(p.algo) + ".csv"));

    csv::Document dist_out;
    utils::AppendColumns(dist_out, {"avg_distance"}, std::make_tuple(distances));
    dist_out.Save(p.out_dir / (ex09::tag_from(p.algo) + "_stats.csv"));

    return 0;
}
