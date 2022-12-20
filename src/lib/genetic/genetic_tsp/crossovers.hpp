//
// Created by Davide Nicoli on 23/11/22.
//

#ifndef ESERCIZI_LSN_GENETIC_TSP_CROSSOVERS_HPP
#define ESERCIZI_LSN_GENETIC_TSP_CROSSOVERS_HPP

#include <cstddef>
#include <iterator>
#include <random>

#include "../genetic_utils.hpp"

template<std::ptrdiff_t N_CITIES>
class [[maybe_unused]] ExerciseCrossover {
    static constexpr const std::ptrdiff_t I_SIZE = N_CITIES - 1;

public:
    /**
     * Performs the crossover of two individuals as described by the exercise instructions.
     * @param parent_1 First parent.
     * @param parent_2 Second parent.
     * @param child_1 First child.
     * @param child_2 Second child.
     * @param rng Uniform Random Bit Generator, as specified by the C++ standard.
     */
    template<typename Individual, class URBG>
    void crossover(const Individual &parent_1, const Individual &parent_2, Individual &child_1,
                   Individual &child_2, URBG &rng) {
        using std::next;
        const auto cut = m_cutdist(rng);
        // Orders the elements after "cut" of parent_1 as they appear in parent_2 and vice-versa.
        cut_and_mix(parent_1.cbegin(), parent_2.cbegin(), child_1.begin(), child_2.begin(), I_SIZE,
                    cut, cut, I_SIZE - cut);
    }

private:
    std::uniform_int_distribution<std::ptrdiff_t> m_cutdist{0, N_CITIES - 2};
};

template<std::ptrdiff_t N_CITIES>
class [[maybe_unused]] ExerciseModCrossover {
    using uip = std::uniform_int_distribution<std::ptrdiff_t>::param_type;
    static constexpr const std::ptrdiff_t I_SIZE = N_CITIES - 1;
    static constexpr const std::ptrdiff_t MIN_SLICE_LENGTH = 2;

public:
    /**
     * Performs the crossover of two individuals.
     * @param parent_1 First parent.
     * @param parent_2 Second parent.
     * @param child_1 First child.
     * @param child_2 Second child.
     * @param rng Uniform Random Bit Generator, as specified by the C++ standard.
     */
    template<typename Individual, class URBG>
    void crossover(const Individual &parent_1, const Individual &parent_2, Individual &child_1,
                   Individual &child_2, URBG &rng) {
        using std::next;
        const auto start_slice1 = m_cutdist(rng, m_firstp);
        const auto start_slice2 = m_cutdist(rng, m_firstp);
        const uip length_params{MIN_SLICE_LENGTH,
                                std::min({I_SIZE - start_slice1, I_SIZE - start_slice2})};
        // Orders the "cut" elements of parent_1 as they appear in parent_2 and vice-versa.
        cut_and_mix(parent_1.cbegin(), parent_2.cbegin(), child_1.begin(), child_2.begin(), I_SIZE,
                    start_slice1, start_slice2, m_cutdist(rng, length_params));
    }

private:
    static_assert(N_CITIES >= 1 + MIN_SLICE_LENGTH);
    std::uniform_int_distribution<std::ptrdiff_t> m_cutdist{};
    const uip m_firstp{0, I_SIZE - MIN_SLICE_LENGTH};
};

template<std::ptrdiff_t N_CITIES>
class [[maybe_unused]] MyCrossover2 {
    using uip = std::uniform_int_distribution<std::ptrdiff_t>::param_type;
    static constexpr const std::ptrdiff_t I_SIZE = N_CITIES - 1;
    static constexpr const std::ptrdiff_t MAX_SLICE_LENGTH = I_SIZE - 1;
    static constexpr const std::ptrdiff_t MIN_SLICE_LENGTH = std::max(0L, MAX_SLICE_LENGTH - 6);

public:
    /**
     * Performs the crossover of two individuals copying slices of one parent to the other.
     * The elements of the transfered slice are removed from their original position to produce valid
     * individuals.
     * @param parent_1 First parent.
     * @param parent_2 Second parent.
     * @param child_1 First child.
     * @param child_2 Second child.
     * @param rng Uniform Random Bit Generator, as specified by the C++ standard.
     */
    template<typename Individual, class URBG>
    void crossover(const Individual &parent_1, const Individual &parent_2, Individual &child_1,
                   Individual &child_2, URBG &rng) {
        // Parameters to build a distribution which samples candidate crossover cuts.
        // These must range
        // Sampling the crossover slice starts for both parents
        const auto cut_start1 = m_cut_dist(rng, m_firstp);
        const auto cut_start2 = m_cut_dist(rng, m_firstp);
        // Parameters to sample a valid slice size
        const uip length_params{
                MIN_SLICE_LENGTH,
                std::min({I_SIZE - cut_start1, I_SIZE - cut_start2, MAX_SLICE_LENGTH})};
        // Sampling the slice size
        const auto slices_length = m_cut_dist(rng, length_params);
        utils::exchange_slices(parent_1.cbegin(), parent_2.cbegin(), child_1.begin(),
                               child_2.begin(), I_SIZE, cut_start1, cut_start2, slices_length);
    }

private:
    static_assert(I_SIZE > MIN_SLICE_LENGTH);
    std::uniform_int_distribution<std::ptrdiff_t> m_cut_dist{};
    const uip m_firstp{0, I_SIZE - MIN_SLICE_LENGTH};
};


struct [[maybe_unused]] CloneCrossover {
    template<typename Individual, class URBG>
    inline void crossover(const Individual &parent_1, const Individual &parent_2,
                          Individual &child_1, Individual &child_2, URBG & /*rng*/) {
        child_1 = parent_2;
        child_2 = parent_1;
    }
};

template<std::ptrdiff_t N_CITIES>
class [[maybe_unused]] Fusion {
public:
    explicit Fusion(double p) : m_coin(p) {}
    template<typename Individual, class URBG>
    void crossover(const Individual &parent_1, const Individual &parent_2, Individual &child_1,
                   Individual &child_2, URBG &rng) {
        if (m_coin(rng)) m_my2.template crossover(parent_1, parent_2, child_1, child_2, rng);
        else
            m_exmod.template crossover(parent_1, parent_2, child_1, child_2, rng);
    }

private:
    MyCrossover2<N_CITIES> m_my2{};
    ExerciseModCrossover<N_CITIES> m_exmod{};
    std::bernoulli_distribution m_coin;
};


#endif//ESERCIZI_LSN_GENETIC_TSP_CROSSOVERS_HPP
