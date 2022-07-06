#ifndef UTILS_HPP
#define UTILS_HPP

#include <algorithm>
#include <functional>
#include <iterator>
#include <numeric>

using std::next;

namespace utils {

    template<typename It>
    inline constexpr auto snext(It it, size_t N) {
        using diff_t = typename std::iterator_traits<It>::difference_type;
        return std::next(it, static_cast<diff_t>(N));
    }

    template<typename InputIt, typename OutputIt>
    constexpr OutputIt partial_average(InputIt first, InputIt last, OutputIt d_first) {
        if (first == last) { return d_first; }

        using Value = typename std::iterator_traits<InputIt>::value_type;
        Value sum = *first;
        *d_first = sum;
        size_t i = 1;
        while (++first != last) {
            sum += *first;
            *++d_first = sum / Value(++i);
        }
        return ++d_first;
    }

    template<typename InputIt, typename OutputIt, typename BinaryOperation>
    constexpr OutputIt blocks(InputIt first, InputIt last, OutputIt d_first, size_t block_size,
                              BinaryOperation binary_op) {
        if (first == last) { return d_first; }
        assert(std::distance(first, last) % block_size == 0);
        while (first != last) {
            *d_first++ = std::accumulate(first, next(first, block_size), binary_op);
            first = next(first, block_size);
        }
        return d_first;
    }

    template<typename InputIt, typename OutputIt>
    constexpr OutputIt blocks(InputIt first, InputIt last, OutputIt d_first, size_t block_size) {
        if (first == last) { return d_first; }
        assert(std::distance(first, last) % block_size == 0);
        while (first != last) {
            *d_first++ = std::reduce(first, next(first, block_size)) / block_size;
            first = next(first, block_size);
        }
        return d_first;
    }

    template<typename InputIt, typename OutputIt, typename Compare>
    constexpr void argsort(InputIt first, InputIt last, OutputIt d_first, Compare comp) {
        const auto N = std::distance(first, last);
        const auto d_last = next(d_first, N);
        std::iota(d_first, d_last, 0);
        std::stable_sort(d_first, d_last, [&](const auto left, const auto right) -> bool {
            // sort indices according to corresponding array element
            return comp(*next(first, left), *next(first, right));
        });
    }

    template<typename InputIt, typename OutputIt>
    constexpr void argsort(InputIt first, InputIt last, OutputIt d_first) {
        argsort(first, last, d_first, std::less<>());
    }

    template<typename Vector>
    constexpr typename Vector::value_type norm2(const Vector &v) {
        return std::transform_reduce(v.cbegin(), v.cend(), v.cbegin(),
                                     typename Vector::value_type(0));
    }

    template<typename Vector>
    constexpr typename Vector::value_type distance2(const Vector &from, const Vector &to) {
        return std::transform_reduce(
                from.cbegin(), from.cend(), to.cbegin(), 0, std::plus<>(),
                [](const auto x, const auto y) { return std::pow((y - x), 2); });
    }

    // template <typename InputIt>
    // void sort
}// namespace utils

#endif
