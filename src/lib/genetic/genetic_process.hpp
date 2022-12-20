//
// Created by Davide Nicoli on 19/05/22.
//

#ifndef GENETIC_TSP_GENETIC_PROCESS_HPP
#define GENETIC_TSP_GENETIC_PROCESS_HPP

#include "config.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <iostream>
#include <random>
#include <vector>

#if USE_MPI
#include <mpi.h>
#endif

#include <indicators/dynamic_progress.hpp>
#include <indicators/progress_bar.hpp>
#include <rapidcsv.h>

#include "genetic_utils.hpp"
#include "utils.hpp"

using namespace utils;
namespace csv = rapidcsv;

namespace genetic {
    template<size_t N_CITIES>
    void load_coordinates(const fs::path &path,
                          std::array<std::array<double, 2>, N_CITIES> &coordinates, bool header) {
        const int header_param = header ? 0 : -1;
        csv::Document table(path, rapidcsv::LabelParams(header_param, -1));
        std::vector<double> buffer;
        for (size_t i = 0; i < N_CITIES; i++) {
            buffer = table.GetRow<double>(i);
            std::copy(buffer.cbegin(), buffer.cend(), coordinates[i].begin());
        }
    }

    /**
     * Class describing a genetic process. It accepts a genetic algorithm as a template parameter
     * which describe what an individual is and how mutation, crossover, and evaluation are performed.
     * @tparam GA Genetic Algorithm. Must implement certain methods (see genetic_tspalgorithms/tsp_ga.hpp).
     */
    template<class GA>
    class Process {

    public:
        typedef typename GA::Individual Individual;
        explicit Process(GA &&ga) : m_ga(std::forward<GA>(ga)) {}

        /**
         * Generates a random population
         * @param first_individual First element of the population container.
         * @param N Number of individuals in the population.
         * @param rng Uniform random bit generator, as defined by the C++ standard.
         */
        template<typename PopulationIt, class URBG>
        inline constexpr void generate(PopulationIt first_individual, size_t N, URBG &rng) {
            m_ga.generate(first_individual, N, rng);
        }

        /**
         * Computes the fitness of each individual.
         * @param first_individual Iterator to the first element of the population.
         * @param N Number of individuals in the population.
         * @param first_evaluation Iterator to the first element of the evaluations container.
         */
        template<typename PopulationIt, typename EvaluationsIt>
        inline constexpr void evaluate(PopulationIt first_individual, size_t N,
                                       EvaluationsIt first_evaluation) {
            std::transform(first_individual, utils::snext(first_individual, N), first_evaluation,
                           [&](const auto &i) { return m_ga.fitness(i); });
        }

        /**
         * Selects parents from a population using their fitness and stores them in a container.
         * @param first_individual Iterator to the first individual.
         * @param N Number of individuals and parents.
         * @param first_parent Iterator to the first element of the parents' container.
         * @param first_evaluation Iterator to the first element of the fitnesses' container.
         * @param rng Uniform random bit generator, as defined by the C++ standard.
         */
        template<typename PopulationIt, typename ParentIt, typename EvaluationsIt, class URBG>
        inline constexpr void select_parents(PopulationIt first_individual, size_t N,
                                             ParentIt first_parent, EvaluationsIt first_evaluation,
                                             URBG &rng) {
            m_ga.select_parents(first_individual, N, first_parent, first_evaluation, rng);
        }

        /**
         * Performs the crossover operation on a list of parents and stores the results in another container.
         * @param first_parent Iterator to the first parent.
         * @param N Number of parents and children.
         * @param first_child Iterator to the first element of the children' container.
         * @param rng Uniform random bit generator, as defined by the C++ standard.
         */
        template<typename ParentIt, typename ChildIt, class URBG>
        inline constexpr void crossover(ParentIt first_parent, size_t N, ChildIt first_child,
                                        URBG &rng) {
            for (size_t i = 0; i < N; i += 2) {
                m_ga.crossover(*first_parent++, *first_parent++, *first_child++, *first_child++,
                               rng);
            }
        }

        /**
         * Performs the mutation operation at random.
         * @param first_individual Iterator to the first individual.
         * @param N Number of individuals and parents.
         * @param mutation_probability Each individual's probability of undergoing a mutation.
         * @param rng Uniform random bit generator, as defined by the C++ standard.
         */
        template<typename PopulationIt, class URBG>
        inline constexpr void mutate(PopulationIt first_individual, size_t N,
                                     double mutation_probability, URBG &rng) {
            std::for_each(first_individual, snext(first_individual, N), [&](auto &individual) {
                if (m_mutprob(rng) < mutation_probability) m_ga.mutate(individual, rng);
            });
        }

#if USE_MPI
        /**
         * Performs the standard genetic process evolution loop: generation -> evaluation -> parent selection -> crossover -> mutation -> evaluation -> ...
         * The evolution is split in blocks of generations for MPI's sake: at the end of each block
         * the population is synchronized choosing the (population_size / n_processes) best scoring
         * individuals. Note: in the end only the first process gathers the best individuals.
         * @param first_individual Iterator to the first individual.
         * @param population_size Population size.
         * @param first_evaluation Output iterator to the first element of the fitnesses' container.
         * @param migration_length Number of iterations for each block.
         * @param n_migrations Number of blocks (or number of synchronizations).
         * @param mutation_probability Mutation probability.
         * @param rng Uniform random bit generator, as defined by the C++ standard.
         */
        template<typename PopulationIt, typename EvaluationsIt, class URBG>
        void mpi_run(PopulationIt first_individual, size_t population_size,
                     EvaluationsIt first_evaluation, size_t migration_length, size_t n_migrations,
                     double mutation_probability, indicators::ProgressBar &pbar, URBG &rng) {
            using namespace indicators;

            if (n_migrations == 0) return;

            if (population_size % 2 == 1) {
                throw std::runtime_error(
                        "Population must be composed of an even number of individuals");
            }

            // MPI initialization
            int mpi_id;
            int n_procs;
            MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);
            MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
            // The popoulation size must be divisible by the number of processes, to avoid complications.
            if (population_size % size_t(n_procs) != 0UL) {
                throw std::runtime_error(
                        "Population size should be a multiple of the number of processes.\n"
                        "population_size: " +
                        std::to_string(population_size) + "\tn_procs: " + std::to_string(n_procs));
            }
            // The number of individuals which will be shared with the other processes
            const auto n_per_process = signed(population_size) / n_procs;
            // Vector where the global rank of each individual will be stored
            std::vector<size_t> ranks(population_size);
            std::vector<Individual> parents_buffer(population_size);

            // Populating and evaluating the first generation.
            generate(first_individual, population_size, rng);
            evaluate(first_individual, population_size, first_evaluation);

            for (auto i = 0U; i < n_migrations; i++) {
                // Performing the select_parents -> crossover -> mutate -> fitness loop for "migration_length" times.
                migration_loop(first_individual, population_size, parents_buffer.begin(),
                               first_evaluation, migration_length, mutation_probability, rng);
                // Gathering the best "n_per_process" individuals from each process into parents_buffer
                gather_best_individuals(first_individual, population_size, parents_buffer.data(),
                                        first_evaluation, ranks.begin(), n_per_process,
                                        (i + 1 < n_migrations));
                std::copy(parents_buffer.cbegin(), parents_buffer.cend(), first_individual);
                // Could also be gathered for better performance
                evaluate(first_individual, population_size, first_evaluation);
                if (mpi_id == 0) {
                    // Prints the progress bar and shows the best fitness.
                    update_best(first_individual, population_size, first_evaluation);
                    pbar.set_option(
                            option::PostfixText{"best fitness: " + std::to_string(m_best_fitness)});
                    pbar.tick();
                }
            }
        }
#endif
        /**
         * Performs the standard genetic process evolution loop: generation -> evaluation -> parent selection -> crossover -> mutation -> evaluation -> ...
         * The evolution is split in blocks of generations for MPI's sake: at the end of each block
         * the population is synchronized choosing the (population_size / n_processes) best scoring
         * individuals.
         * @param first_individual Iterator to the first individual.
         * @param population_size Population size.
         * @param first_evaluation Output iterator to the first element of the fitnesses' container.
         * @param n_iterations Number of iterations.
         * @param mutation_probability Mutation probability.
         * @param rng Uniform random bit generator, as defined by the C++ standard.
         */
        template<bool compute_statistic, typename PopulationIt, typename EvaluationsIt,
                 typename StatisticsIt, class URBG>
        void run(PopulationIt first_individual, size_t population_size,
                 EvaluationsIt first_evaluation, size_t n_iterations, double mutation_probability,
                 StatisticsIt first_statistics, indicators::ProgressBar &pbar, URBG &rng) {
            using namespace indicators;

            if (population_size % 2 == 1) {
                throw std::runtime_error(
                        "Population must be composed of an even number of individuals");
            }

            // Buffer where chosen parents are stored in between generations
            std::vector<Individual> parents_buffer(population_size);
            // Populating and evaluating the first generation
            generate(first_individual, population_size, rng);
            evaluate(first_individual, population_size, first_evaluation);

            for (auto i = 0U; i < n_iterations; i++) {
                // Selecting the best scoring individuals
                select_parents(first_individual, population_size, parents_buffer.begin(),
                               first_evaluation, rng);
                if constexpr (compute_statistic) {
                    // Computing the statistics. The funcion is destructive of the fitnesses: it is placed here
                    // as they are not needed anymore and the next line overwrites them anyway
                    *first_statistics++ = unsafe_fitness_statistics(
                            first_evaluation, population_size, [&](const auto fitness) {
                                return m_ga.statistic_from_fitness(fitness);
                            });
                }
                // Performing crossover -> mutation -> evaluation
                cross_mut_eval(first_individual, population_size, parents_buffer.cbegin(),
                               first_evaluation, mutation_probability, rng);
                update_best(first_individual, population_size, first_evaluation);
                // Printing the progress bar and showing the best fitness
                pbar.set_option(
                        option::PostfixText{"best fitness: " + std::to_string(m_best_fitness)});
                pbar.tick();
            }
        }

        auto get_best() const { return std::make_tuple(m_best_overall, m_best_fitness); }

    private:
        GA m_ga;
        std::uniform_real_distribution<double> m_mutprob{};
        Individual m_best_overall;
        double m_best_fitness{0};

#if USE_MPI
        /**
         * Gathers the best individuals found by each MPI process.
         * @param first_individual Iterator to the first individual of the running process.
         * @param population_size Population size.
         * @param first_buffer Iterator to the first element of a buffer of population_size, where the best "n_per_process" individuals of each process will be syncronized.
         * @param first_evaluation Output iterator to the first element of the fitnesses' container.
         * @param first_rank Output iterator to the first element of the ranks' container.
         * @param n_per_process Number of individuals per process.
         * @param copy_to_all_processes Whether to syncronize the buffer of every process or just the first one's.
         */
        template<typename PopulationIt, typename BufferIt, typename EvaluationsIt, typename RanksIt>
        inline void gather_best_individuals(PopulationIt first_individual, size_t population_size,
                                            BufferIt first_buffer, EvaluationsIt first_evaluation,
                                            RanksIt first_rank, int n_per_process,
                                            bool copy_to_all_processes) const {
            // Computing the individual's ranks placing first the elements of greatest fitness
            rank_n(first_evaluation, population_size, first_rank, std::greater<>());
            // Sorting the population depending on that
            order_to_n(first_individual, population_size, first_rank);

            if (copy_to_all_processes)
                // Synchronizing buffers of different processes. This gathers the
                // best ranking individuals for each process, not the best ranking among processes.
                // This way there is more entropy.
                MPI_Allgather(&*first_individual, n_per_process, m_ga.individual_mpi, first_buffer,
                              n_per_process, m_ga.individual_mpi, MPI_COMM_WORLD);
            else
                // As above, but the indivuals are copied only to the first process' buffer for efficiency
                MPI_Gather(&*first_individual, n_per_process, m_ga.individual_mpi, first_buffer,
                           n_per_process, m_ga.individual_mpi, 0, MPI_COMM_WORLD);
        }
#endif

        /**
         * Grouping crossover -> mutation -> evaluation in a single command.
         */
        template<typename PopulationIt, typename ParentIt, typename EvaluationsIt, class RNG>
        inline void cross_mut_eval(PopulationIt first_individual, size_t population_size,
                                   ParentIt first_parent, EvaluationsIt first_evaluation,
                                   double mutation_probability, RNG &rng) {
            crossover(first_parent, population_size, first_individual, rng);
            mutate(first_individual, population_size, mutation_probability, rng);
            evaluate(first_individual, population_size, first_evaluation);
        }

        template<typename PopulationIt, typename BufferIt, typename EvaluationsIt, class RNG>
        inline void migration_loop(PopulationIt first_individual, size_t population_size,
                                   BufferIt first_buffer, EvaluationsIt first_evaluation,
                                   size_t n_iterations, double mutation_probability, RNG &rng) {
            for (size_t i = 0UL; i < n_iterations; i++) {
                select_parents(first_individual, population_size, first_buffer, first_evaluation,
                               rng);
                cross_mut_eval(first_individual, population_size, first_buffer, first_evaluation,
                               mutation_probability, rng);
                update_best(first_individual, population_size, first_evaluation);
            }
        }

        template<typename PopulationIt, typename EvaluationsIt>
        inline void update_best(PopulationIt first_individual, size_t population_size,
                                EvaluationsIt first_evaluation) {
            // Printing the progress bar and showing the best fitness
            const auto best_fitness = std::max_element(
                    first_evaluation, utils::snext(first_evaluation, population_size));
            // Updating the best element
            if (*best_fitness > m_best_fitness) {
                m_best_overall =
                        *std::next(first_individual, std::distance(first_evaluation, best_fitness));
                m_best_fitness = *best_fitness;
            }
        }

        /**
         * WATCH OUT: this function mutates the evaluations and quite specific to the TSP. Should be reworked.
         * @param first_evaluation
         * @param population_size
         * @return
         */
        template<typename EvaluationsIt, class Statistic>
        double unsafe_fitness_statistics(EvaluationsIt first_evaluation, size_t population_size,
                                         Statistic stat_from_fitness) {
            auto median_pointer = utils::snext(first_evaluation, population_size / 2);
            // O(pop_size)
            std::nth_element(first_evaluation, median_pointer,
                             utils::snext(first_evaluation, population_size));
            double stat_sum = 0.0;
            size_t tot_elements = 0;
            // O(pop_size)
            for (size_t i = 0UL; i < population_size; i++) {
                if (*first_evaluation >= *median_pointer) {
                    stat_sum += stat_from_fitness(*first_evaluation);
                    tot_elements++;
                }
                first_evaluation++;
            }
            return stat_sum / double(tot_elements);
        }
    };
}// namespace genetic

#endif// GENETIC_TSP_GENETIC_PROCESS_HPP
