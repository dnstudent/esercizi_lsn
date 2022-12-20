//
// Created by Davide Nicoli on 30/05/22.
//

#ifndef GENETIC_TSP_TSP_GA_HPP
#define GENETIC_TSP_TSP_GA_HPP

#include "config.hpp"

#include <array>
#include <numeric>
#include <random>
#include <utility>
#include <vector>

#if USE_MPI
#include <mpi.h>
#endif

#include "../genetic_utils.hpp"
#include "utils.hpp"

using namespace utils;

/**
 * Methods and structures to cope with the Travelling Salesman Problem (where each city is visited only once) in this framework.
 * @tparam Coordinates Type name for the coordinates points.
 * @tparam N_CITIES Number of cities, defined at compile time for performance reasons.
 */
template<class Crossover, typename Coordinates, size_t N_CITIES>
class TSP {
    typedef uint16_t city_index;
    // Preventing numeric and stack overflow at compile time
    static_assert((N_CITIES <= std::min(1000UL, size_t(std::numeric_limits<city_index>::max()))));
    static constexpr const size_t I_SIZE{N_CITIES - 1};

public:
    /**
     * Each individual is an ordered collection of city indices. The first one is fixed at 0.
     */
    typedef std::array<city_index, I_SIZE> Individual;
    typedef double FitnessMeasure;

    TSP(const std::array<Coordinates, N_CITIES> &city_coordinates, Crossover &&cross_fn)
        : m_city_coordinates(city_coordinates), m_crossover(std::forward<Crossover>(cross_fn)) {
#if USE_MPI
        MPI_Type_contiguous(I_SIZE, MPI_UNSIGNED_SHORT, &individual_mpi);
        MPI_Type_commit(&individual_mpi);
#endif
    }

    /**
     * Generate a population of Individuals.
     * @param first_individual Output iterator to the first individual.
     * @param N Population size.
     * @param rng Uniform Random Bit Generator as defined by the C++ standard.
     */
    template<typename PopulationIt, class URBG>
    void generate(PopulationIt first_individual, size_t N, URBG &rng) {
        Individual base;

        std::iota(base.begin(), base.end(), 1);
        std::for_each(first_individual, utils::snext(first_individual, N), [&](Individual &i) {
            std::copy(base.cbegin(), base.cend(), i.begin());
            std::shuffle(i.begin(), i.end(), rng);
        });
    }

    /**
     * Computes an individual's fitness.
     * @param individual Individual.
     * @return Fitness.
     */
    FitnessMeasure fitness(const Individual &individual) {
        // The first city is fixed
        auto total_distance = std::sqrt(
                utils::distance2(m_city_coordinates[0], m_city_coordinates[*individual.cbegin()]));
        for (auto i = std::next(individual.cbegin()); i < individual.cend(); i++) {
            total_distance += std::sqrt(
                    utils::distance2(m_city_coordinates[*i], m_city_coordinates[*std::prev(i)]));
        }
        return static_cast<FitnessMeasure>(1.0) / static_cast<FitnessMeasure>(total_distance);
    }

    auto statistic_from_fitness(const FitnessMeasure fitness) {
        return static_cast<FitnessMeasure>(1.0) / fitness;
    }

    /**
     * Selects the individuals most suited to be parents from a population.
     * @param first_individual Iterator to the first individual of a population.
     * @param N Population size.
     * @param first_new_individual Output iterator to the parents' container.
     * @param first_evaluation Input iterator to the fitness' container.
     * @param rng Uniform Random Bit Generator, as specified by the C++ standard.
     */
    template<typename PopulationIt, typename EvaluationsIt, typename OutPopIt, class URBG>
    static void select_parents(PopulationIt first_individual, size_t N,
                               OutPopIt first_new_individual, EvaluationsIt first_evaluation,
                               URBG &rng) {
        std::discrete_distribution<int64_t> parents_distribution(first_evaluation,
                                                                 utils::snext(first_evaluation, N));
        std::generate_n(first_new_individual, N,
                        [&]() { return *std::next(first_individual, parents_distribution(rng)); });
    }


    /**
     * Performs the crossover of two individuals. Proxy for different methods.
     * @tparam dummy Whether to use a crossover algorithm or
     */
    template<class URBG>
    inline void crossover(const Individual &parent_1, const Individual &parent_2,
                          Individual &child_1, Individual &child_2, URBG &rng) {
        return m_crossover.crossover(parent_1, parent_2, child_1, child_2, rng);
    }

    /**
     * Performs a randomly chosen mutation on an individual.
     * @param individual The individual.
     * @param rng Uniform Random Bit Generator, as specified by the C++ standard.
     */
    template<class URBG>
    void mutate(Individual &individual, URBG &rng) {
        const auto roll = m_mutation_distribution(rng);
        if (roll == 0) {
            _mutate_reflect(individual, rng);
        } else if (roll == 1) {
            _mutate_swap_ranges(individual, rng);
        } else if (roll == 2) {
            _mutate_shift(individual, rng);
        }
    }

#if USE_MPI
    MPI_Datatype individual_mpi{};
#endif

private:
    const std::array<Coordinates, N_CITIES> m_city_coordinates;
    Crossover m_crossover;
    using uip = std::uniform_int_distribution<std::ptrdiff_t>::param_type;
    // Parameters to build a distribution which samples candidate starting cuts.
    // Distribution parameters: prevent picking the last city, otherwise the cut produce no effect.
    const uip m_startp{0, I_SIZE - 1};
    // Dsitribution paramters to sample candidate ending cuts.
    // Prevent picking the first city, otherwise the cut produce no effect.
    const uip m_endp{1, I_SIZE};
    // Distribution object which will be used to propose candidate cut points.
    std::uniform_int_distribution<std::ptrdiff_t> m_cut_dist{};
    // Mutation dice.
    std::uniform_int_distribution<unsigned short> m_mutation_distribution{0, 2};


    /**
     * Randomly chooses a slice of an individual and reflects its elements.
     */
    template<class URBG>
    void _mutate_reflect(Individual &individual, URBG &rng) {
        // Samples two indices to elements inside the individual
        const auto i1 = m_cut_dist(rng, m_startp);
        const auto i2 = m_cut_dist(rng, uip{i1 + 1, I_SIZE});
        //        if (i1 > i2) { std::swap(i1, i2); }
        std::reverse(std::next(individual.begin(), i1), std::next(individual.begin(), i2));
    }

    /**
     * Randomly chooses two equally long slices of an individual and swaps them.
     */
    template<class URBG>
    void _mutate_swap_ranges(Individual &individual, URBG &rng) {
        using std::next;
        // Samples four possible cut points inside the individual
        std::array<int, 4> cuts{};
        std::generate(cuts.begin(), cuts.begin() + 2, [&]() { return m_cut_dist(rng, m_startp); });
        std::generate(cuts.begin() + 2, cuts.end(), [&]() { return m_cut_dist(rng, m_endp); });
        // sorts them so that, taken in sequence, they determine two valid and non-overlapping slices
        std::sort(cuts.begin(), cuts.end());
        const auto first = individual.begin();
        // taking the shortest slice length
        const auto length = std::min(cuts[1] - cuts[0], cuts[3] - cuts[2]);
        // inverting the slices
        std::swap_ranges(next(first, cuts[0]), next(first, cuts[0] + length), next(first, cuts[2]));
    }

    /**
     * Shifts a random slice by a random number of postions.
     */
    template<class URBG>
    void _mutate_shift(Individual &individual, URBG &rng) {
        using std::next;
        // Generating the cut points
        static const uip cut_params{0, I_SIZE};
        std::array<std::ptrdiff_t, 3> cuts{};
        std::generate(cuts.begin(), cuts.end(), [&]() { return m_cut_dist(rng, cut_params); });
        // Sorting them so that they determine two contiguous slices
        std::sort(cuts.begin(), cuts.end());
        // Performing the shift
        const auto first = individual.begin();
        std::rotate(next(first, cuts[0]), next(first, cuts[1]), next(first, cuts[2]));
    }
};

#endif// GENETIC_TSP_TSP_GA_HPP
