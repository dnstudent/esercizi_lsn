#ifndef ESTIMATORS_HPP
#define ESTIMATORS_HPP

#include <array>
#include <cstddef>
#include <tuple>
#include <utility>
#include <vector>

#include "../utils.hpp"

namespace estimators {
    /**
     * Computes and stores variables estimates from transform dataset
     * @param first Iterator to first data entry.
     * @param last Iterator past last data entry.
     * @param estimators Tuple of estimators.
     * @param outs Iterators/pointers to results storage. For each estimator there must be two output pointers (the first one for the estimate, the second one for the uncertainty).
     * @return
     */
    template<size_t I = 0, typename SampleIt, class Estimators, typename... OutputIt>
    constexpr void store_estimates(SampleIt first, SampleIt last, Estimators &estimators,
                                   OutputIt... outs) {
        static_assert(2 * std::tuple_size_v<Estimators> == sizeof...(outs));
        auto outs_tuple = std::tie(outs...);
        if constexpr (I == std::tuple_size_v<Estimators>) {
            return;
        } else {
            std::tie(*std::get<2 * I>(outs_tuple), *std::get<2 * I + 1>(outs_tuple)) =
                    std::get<I>(estimators)(first, last);
            store_estimates<I + 1>(first, last, estimators, outs...);
        }
    }

}// namespace estimators

#endif