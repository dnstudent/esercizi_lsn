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
    template<typename value>
    class Variance {
    public:
        explicit Variance(const value mean) : m_mean{mean} {}

        template<typename It>
        constexpr std::pair<value, value> operator()(It first, It last) const {
            return Mean<value>()(first, last,
                                 [&](const auto x) { return (x - m_mean) * (x - m_mean); });
        }

    private:
        const value m_mean;
    };

    /**
     * Estimates a variable variance using sample variance with a given mean. Computation is cached so that results constitute a progressive estimate. Provides the statistical uncertainty.
     * @tparam value The numeric field to use.
     */
    template<typename value>
    class StatefulVariance {
    public:
        explicit StatefulVariance(const value mean) : m_mean{mean} {}

        /**
         * Stateful estimation of the variance with uncertainty from a variable sample.
         * @param first Iterator to the first element of the sample (sample.begin()).
         * @param last Iterator past the last element of the sample (sample.end()).
         * @return Pair (estimate, uncertainty)
         */
        template<typename It>
        constexpr auto operator()(It first, It last) {
            const auto block_size = std::distance(first, last);
            m_current_block++;
            // TODO: if (block_size == 0)...
            m_transformed.resize(block_size);
            std::transform(first, last, m_transformed.begin(), [&](const auto x) {
                const auto u = x - m_mean;
                return u * u;
            });
            const auto block_mean =
                    std::reduce(m_transformed.cbegin(), m_transformed.cend(), value(0)) /
                    value(block_size);
            m_running_sum += block_mean;
            const auto mean_estimate = m_running_sum / value(m_current_block);
            m_running_sum2 += block_mean * block_mean;

            if (m_current_block == 1) { return std::make_pair(mean_estimate, value(0)); }

            const auto square_mean_estimate = m_running_sum2 / value(m_current_block);
            const value variance_estimate = square_mean_estimate - mean_estimate * mean_estimate;
            return std::make_pair(mean_estimate,
                                  std::sqrt(variance_estimate / value(m_current_block - 1)));
        }

    private:
        const value m_mean;
        size_t m_current_block{0};
        value m_running_sum{0};
        value m_running_sum2{0};
        std::vector<value> m_transformed{};
    };

}// namespace estimators

#endif//ESERCIZI_LSN_VARIANCE_HPP
