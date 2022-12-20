//
// Created by Davide Nicoli on 06/07/22.
//

#ifndef ESERCIZI_LSN_CHI2_HPP
#define ESERCIZI_LSN_CHI2_HPP

#include <functional>
#include <numeric>

namespace estimators {
    template<typename value>
    class UniformChi2 {
    public:
        explicit UniformChi2(value expected_value) : m_expected_value{expected_value} {}

        typedef value Output;

        /**
         * Computes the X^2 statistic for transform sample
         * @param first The first element of the sample.
         * @param last The last element of the sample.
         * @return X^2
         */
        template<typename It>
        constexpr Output operator()(It first, It last) {
            return std::transform_reduce(first, last, value(0), std::plus<>(), [&](const auto Oi) {
                auto num = value(Oi) - m_expected_value;
                return num * num / m_expected_value;
            });
        }

    private:
        value m_expected_value;
    };
}// namespace estimators

#endif//ESERCIZI_LSN_CHI2_HPP
