//
// Created by Davide Nicoli on 04/07/22.
//

#ifndef ESERCIZI_LSN_VARIANCE_HPP
#define ESERCIZI_LSN_VARIANCE_HPP

#include <array>
#include <cmath>
#include <functional>
#include <iterator>
#include <numeric>
#include <utility>

#include "../utils.hpp"
#include "mean.hpp"

namespace estimators {
    /**
     * Estimates a variable's variance as the average of rejects with statistical uncertainty.
     * @tparam value The numeric field to use.
     */
    template<typename value>
    class ProgVariance {
    public:
        typedef std::tuple<value, value> Output;
        /**
         * Progressive estimation of a variable's variance with uncertainty.
         * @param first Iterator to the first element of the sample (sample.begin()).
         * @param last Iterator past the last element of the sample (sample.end()).
         * @return Pair (estimate, uncertainty)
         */
        template<typename It>
        constexpr Output operator()(It first, It last) {
            const auto block_size = std::distance(first, last);
            m_current_block++;
            // TODO: if (block_size == 0)...
            m_rejects2.resize(size_t(block_size));
            const auto block_mean = value(std::reduce(first, last)) / value(block_size);
            m_running_mean_sum += block_mean;
            // estimate of the mean up to the current block.
            const auto mean_estimate = m_running_mean_sum / value(m_current_block);
            std::transform(first, last, m_rejects2.begin(), [&](const auto x) {
                const auto s = x - mean_estimate;
                return s * s;
            });
            return m_mean_estimator(m_rejects2.cbegin(), m_rejects2.cend());
        }

    private:
        size_t m_current_block{0};
        ProgAvg<value> m_mean_estimator{};
        value m_running_mean_sum{0};
        std::vector<value> m_rejects2{};
    };

}// namespace estimators

#endif//ESERCIZI_LSN_VARIANCE_HPP
