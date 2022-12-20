//
// Created by Davide Nicoli on 21/03/22.
//

#ifndef ESERCIZI_LSN_DISCRETE_HPP
#define ESERCIZI_LSN_DISCRETE_HPP

#include <functional>
#include <iterator>
#include <map>
#include <random>
#include <vector>

#include "uniform_int.hpp"
#include "utils.hpp"

namespace distributions {
    /**
     * A discrete uniform distribution returning sampled items.
     * @tparam ItemType
     */
    template<class ItemType>
    class UniformDiscrete {
    public:
        typedef const ItemType &result_type;

        template<typename ItemsIt>
        UniformDiscrete(ItemsIt first, ItemsIt last)
            : m_d(0, size_t(std::distance(first, last) - 1)), m_items(first, last) {
            assert(std::distance(first, last) >= 1);
        }
        UniformDiscrete(std::initializer_list<ItemType> items)
            : UniformDiscrete(items.begin(), items.end()) {}

        /**
         * Samples from the uniform discrete distribution
         * @param rng The random engine used for sampling.
         * @return The sampled item.
         */
        template<class URBG>
        result_type operator()(URBG &rng) {
            return m_items[m_d(rng)];
        }

    private:
        uniform_int<size_t> m_d;
        const std::vector<ItemType> m_items;
    };

}// namespace distributions

#endif// ESERCIZI_LSN_DISCRETE_HPP
