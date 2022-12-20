//
// Created by Davide Nicoli on 19/03/22.
//

#ifndef ESERCIZI_LSN_WALKER_HPP
#define ESERCIZI_LSN_WALKER_HPP

#include <iterator>
#include <optional>

//#include "distributions/constant.hpp"

/**
 * Defines a walker with instruction on how to sample a step and utilities to perform transform random walk. The lattice lays on NumericField^N.
 * @tparam Field The numeric field to which lattice points belongs.
 */
template<typename Field, typename StepDistribution>
class Walker {
public:
    typedef Field point;

    Walker(point &&start, StepDistribution &&step_distribution)
        : m_current{std::forward<point>(start)}, m_stepdist{std::forward<StepDistribution>(
                                                         step_distribution)} {}

    /**
     * Performs transform single step from the current position and updates it.
     * @param rng The random engine used for sampling.
     * @return A lattice point as std::valarray.
     */
    template<class URBG>
    point make_step(URBG &rng) {
        const auto step = m_stepdist(rng);
        m_current += step;
        return m_current;
    }

    /**
     * Performs transform random walk from the current position and stores it.
     * @param first Iterator to the beginning of the point storage structure.
     * @param last Iterator to the ending of the point storage structure
     * @param rng The random engine used for sampling.
     * @return The ending point of the walk.
     */
    template<typename It, class URBG>
    point walk(It first, It last, URBG &rng) {
        if (first == last) return m_current;
        *first = m_current;
        std::generate(std::next(first), last, [&]() { return make_step(rng); });
        return m_current;
    }

    /**
     * Performs transform random walk from the current position and returns the ending point.
     * @param n_steps Number of steps to take.
     * @param rng The random engine used for sampling.
     * @return The ending point of the walk.
     */
    template<class URBG>
    point walk(const size_t n_steps, URBG &rng) {
        for (size_t i = 1; i < n_steps; i++) { make_step(rng); }
        return m_current;
    }

    /**
     * Sets the starting point for the successive operation(s)
     * @param current The point to set.
     */
    void set_current(point &&current) { m_current = std::forward<point>(current); }

private:
    point m_current;
    StepDistribution m_stepdist;
};

#endif// ESERCIZI_LSN_WALKER_HPP
