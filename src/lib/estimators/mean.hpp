//
// Created by Davide Nicoli on 04/07/22.
//

#ifndef ESERCIZI_LSN_MEAN_HPP
#define ESERCIZI_LSN_MEAN_HPP

#include <cmath>
#include <functional>
#include <iterator>
#include <numeric>
#include <utility>
#include <vector>

#include "../utils.hpp"

namespace estimators {
    template<typename value>
    struct Mean {
        template<typename It>
        constexpr std::pair<value, value> operator()(It first, It last) const {
            const auto N = std::distance(first, last);
            if (N == 0) { return std::make_pair(0, 0); }
            const auto mean_estimate = std::reduce(first, last, value(0)) / value(N);
            const auto square_mean_estimate =
                    std::transform_reduce(first, last, value(0), std::plus<>(),
                                          [](const auto x) { return x * x; }) /
                    value(N);
            const value variance_estimate = square_mean_estimate - mean_estimate * mean_estimate;
            return std::make_pair(mean_estimate, std::sqrt(variance_estimate / (N - 1)));
        }
        template<typename It, typename Transform>
        constexpr std::pair<value, value> operator()(It first, It last, Transform f) const {
            const auto N = std::distance(first, last);
            if (N == 0) { return std::make_pair(0, 0); }
            std::vector<value> y(N);
            std::transform(first, last, y.begin(), f);
            const auto mean_estimate = std::reduce(y.cbegin(), y.cend(), value(0)) / value(N);
            const auto square_mean_estimate =
                    std::transform_reduce(y.cbegin(), y.cend(), value(0), std::plus<>(),
                                          [](const auto x) { return x * x; }) /
                    value(N);
            const value variance_estimate = square_mean_estimate - mean_estimate * mean_estimate;
            return std::make_pair(mean_estimate, std::sqrt(variance_estimate / (N - 1)));
        }
    };

    /**
     * Estimates a variable mean using sample average. Computation is cached so that results constitute a progressive estimate. Provides the statistical uncertainty.
     * @tparam value The numeric field to use.
     */
    template<typename value>
    class StatefulMean {
    public:
        /**
         * Stateful estimation of the mean with uncertainty from a variable sample.
         * @param first Iterator to the first element of the sample (sample.begin()).
         * @param last Iterator past the last element of the sample (sample.end()).
         * @return Pair (estimate, uncertainty)
         */
        template<typename It>
        constexpr std::pair<value, value> operator()(It first, It last) {
            const auto block_size = std::distance(first, last);
            m_current_block++;
            // TODO: if (block_size == 0)...
            const auto block_mean = std::reduce(first, last, value(0)) / value(block_size);
            m_running_sum += block_mean;
            // <A>
            const auto mean_estimate = m_running_sum / value(m_current_block);
            m_running_sum2 += block_mean * block_mean;

            if (m_current_block == 1) { return {mean_estimate, value(0)}; }

            // <A^2>
            const auto square_mean_estimate = m_running_sum2 / value(m_current_block);
            const value variance_estimate = square_mean_estimate - mean_estimate * mean_estimate;
            return {mean_estimate, std::sqrt(variance_estimate / value(m_current_block - 1))};
        }

        void reset() {
            m_current_block = 0;
            m_running_sum = 0;
            m_running_sum2 = 0;
        }

    protected:
        // Index of the currently processed block
        size_t m_current_block{0};
        // Accumulators for the sample average sum and sum^2
        value m_running_sum{0};
        value m_running_sum2{0};
    };

}// namespace estimators

#endif//ESERCIZI_LSN_MEAN_HPP
