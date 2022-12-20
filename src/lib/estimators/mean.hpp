//
// Created by Davide Nicoli on 04/07/22.
//

#ifndef ESERCIZI_LSN_ESTIMATORS_MEAN_HPP
#define ESERCIZI_LSN_ESTIMATORS_MEAN_HPP

#include <cmath>
#include <functional>
#include <iterator>
#include <numeric>
#include <utility>
#include <vector>

#include "utils.hpp"

namespace estimators {
    template<typename value>
    struct Average {
        typedef std::tuple<value, value> Output;
        template<typename It>
        constexpr Output operator()(It first, It last) const {
            const auto N = std::distance(first, last);
            if (N == 0) { return {0, 0}; }
            const auto mean_estimate = value(std::reduce(first, last)) / value(N);
            const auto square_mean_estimate =
                    value(std::transform_reduce(
                            first, last, typename std::iterator_traits<It>::value_type(0),
                            std::plus<>(), [](const auto x) { return x * x; })) /
                    value(N);
            const value variance_estimate =
                    (square_mean_estimate - mean_estimate * mean_estimate) / value(N - 1);
            return {mean_estimate, std::sqrt(variance_estimate)};
        }
        template<typename It, typename Transform>
        constexpr Output operator()(It first, It last, Transform f) const {
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
     * Estimates a variable mean using progressive sample average. Provides the statistical uncertainty.
     * @tparam value The numeric field to use.
     */
    template<typename value>
    class ProgAvg {
    public:
        typedef std::tuple<value, value> Output;
        /**
         * Progressive estimation of the mean with uncertainty from a variable sample.
         * @param first Iterator to the first element of the sample (sample.begin()).
         * @param last Iterator past the last element of the sample (sample.end()).
         * @return Pair (estimate, uncertainty)
         */
        template<typename It>
        constexpr Output operator()(It first, It last) {
            m_current_block++;
            // TODO: if (block_size == 0)...
            const auto block_avg = utils::average<value>(first, last);
            m_running_sum += block_avg;
            m_running_sum2 += block_avg * block_avg;
            // <A>
            const auto mean_estimate = m_running_sum / static_cast<value>(m_current_block);

            if (m_current_block == 1) { return {mean_estimate, value(0)}; }

            // <A^2>
            const auto square_mean_estimate = m_running_sum2 / static_cast<value>(m_current_block);
            const value estimator_variance =
                    (square_mean_estimate - mean_estimate * mean_estimate) /
                    static_cast<value>(m_current_block - 1);
            return {mean_estimate, std::sqrt(estimator_variance)};
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

    template<typename value>
    class Insta {
    public:
        typedef std::tuple<value, value> Output;
        /**
         * Progressive estimation of the mean with uncertainty from a variable sample.
         * @param first Iterator to the first element of the sample (sample.begin()).
         * @param last Iterator past the last element of the sample (sample.end()).
         * @return Pair (estimate, uncertainty)
         */
        template<typename It>
        constexpr Output operator()(It first, It /**/) {
            return {*first, 0};
        }
    };

    template<typename value>
    class ProgAvg<std::vector<value>> {
        using v = std::vector<value>;

    public:
        explicit ProgAvg(size_t n_bins) : m_estimators(n_bins) {}

        typedef std::tuple<v, v, v> Output;

        template<typename It>
        constexpr Output operator()(It first, It last) {
            Output out{};
            utils::tuple_apply([=](auto &vec) { vec.reserve(m_estimators.size()); }, out);
            for (size_t bin_idx = 0; first != last; bin_idx++, first++) {
                utils::tuple_push_back(m_estimators[bin_idx](first->cbegin(), first->cend()), out);
            }
            return out;
        }

    private:
        std::vector<ProgAvg<value>> m_estimators;
    };

    /**
     * Estimates the mean using sample average. Provides the statistical uncertainty.
     * @tparam value The numeric field to use.
     */
    template<typename value>
    class SampleProgAvg {
    public:
        typedef std::tuple<value, value, value> Output;
        /**
         * Current sample and progressive estimation of the mean with uncertainty.
         * @param first Iterator to the first element of the sample (sample.begin()).
         * @param last Iterator past the last element of the sample (sample.end()).
         * @return Tuple (sample average, progressive estimate, uncertainty)
         */
        template<typename It>
        constexpr Output operator()(It first, It last) {
            m_current_block++;
            // TODO: if (block_size == 0)...
            const auto block_avg = utils::average<value>(first, last);
            m_running_sum += block_avg;
            m_running_sum2 += block_avg * block_avg;
            // <A>
            const auto mean_estimate = m_running_sum / value(m_current_block);

            if (m_current_block == 1) { return {block_avg, mean_estimate, value(0)}; }

            // <A^2>
            const auto square_mean_estimate = m_running_sum2 / value(m_current_block);
            const value estimator_variance =
                    (square_mean_estimate - mean_estimate * mean_estimate) /
                    value(m_current_block - 1);
            return {block_avg, mean_estimate, std::sqrt(estimator_variance)};
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


    /**
     * Here things get reversed for performance reasons
     * @tparam value
     */
    template<typename value>
    class SampleProgAvg<std::vector<value>> {
        using v = std::vector<value>;

    public:
        explicit SampleProgAvg(size_t n_bins) : m_estimators(n_bins) {}

        typedef std::tuple<v, v, v> Output;

        template<typename It>
        constexpr Output operator()(It first, It last) {
            Output out{};
            utils::tuple_apply([=](auto &vec) { vec.reserve(m_estimators.size()); }, out);
            for (size_t bin_idx = 0; first != last; bin_idx++, first++) {
                utils::tuple_push_back(m_estimators[bin_idx](first->cbegin(), first->cend()), out);
            }
            return out;
        }

    private:
        std::vector<SampleProgAvg<value>> m_estimators;
    };

}// namespace estimators

#endif//ESERCIZI_LSN_ESTIMATORS_MEAN_HPP
