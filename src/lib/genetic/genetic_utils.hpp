//
// Created by Davide Nicoli on 25/05/22.
//

#ifndef GENETIC_TSP_UTILS_HPP
#define GENETIC_TSP_UTILS_HPP

#include <algorithm>
#include <cstdlib>
#include <istream>
#include <iterator>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#include "utils.hpp"
namespace utils {
    /**
     * Computes the indices of the elements sorted by compare, as numpy
     * @param first Input iterator to the first element.
     * @param last Input iterator past the last element.
     * @param indices_first Output iterator to the first index.
     * @param compare Operation used to compare elements.
     * @return Iterator past the last index.
     */
    template<typename InputIt, typename OutputIt, typename Compare>
    constexpr auto argsort(InputIt first, InputIt last, OutputIt indices_first, Compare compare) {
        using std::next;
        using diff_t = typename std::iterator_traits<InputIt>::difference_type;
        const auto N = std::distance(first, last);
        const auto indices_last = next(indices_first, N);
        std::iota(indices_first, indices_last, 0);
        return std::sort(
                indices_first, indices_last, [&](const auto left, const auto right) -> bool {
                    // sort indices according to corresponding array element
                    return compare(*next(first, diff_t(left)), *next(first, diff_t(right)));
                });
    }

    /**
     * Same as above but using a number of elements insted of the past-the-end iterator.
     */
    template<typename InputIt, typename OutputIt, typename Compare>
    constexpr inline auto argsort_n(InputIt first, size_t N, OutputIt indices_first,
                                    Compare compare) {
        return argsort(first, snext(first, N), indices_first, compare);
    }

    /**
     * Argsort with a default compare: std::less
     */
    template<typename InputIt, typename OutputIt>
    constexpr inline auto argsort(InputIt first, InputIt last, OutputIt indices_first) {
        return argsort(first, last, indices_first, std::less<>());
    }

    /**
     * As above
     */
    template<typename InputIt, typename OutputIt>
    constexpr inline auto argsort_n(InputIt first, size_t N, OutputIt indices_first) {
        return argsort(first, snext(first, N), indices_first);
    }

    /**
     * Reorders a container placing its elements in the given positions.
     * @param first I/O iterator to the container's first element.
     * @param N Number of elements in the container.
     * @param indices_first Input iterator to the first position. The container must be of length N and its element in the range 0..N-1
     */
    template<typename InputIt, typename OrderIt>
    constexpr inline void order_to_n(InputIt first, size_t N, OrderIt indices_first) {
        const std::vector<typename std::iterator_traits<InputIt>::value_type> elems(
                first, snext(first, N));
        for (auto i = 0UL; i < N; i++) {
            *std::next(first, signed(*snext(indices_first, i))) = elems[i];
        }
    }

    /**
     * Reorders a container placing its elements to the given positions.
     * @param first I/O iterator to the container's first element.
     * @param last I/O iterator to the container's end.
     * @param indices_first Input iterator to the first position. The container must be of length N and its element in the range 0..N-1
     */
    template<typename InputIt, typename OrderIt>
    constexpr inline void order_to(InputIt first, InputIt last, OrderIt indices_first) {
        const auto N = std::distance(first, last);
        return order_to_n(first, unsigned(N), indices_first);
    }

    /**
     * Reorders a container taking elements from the given positions.
     * @param first I/O iterator to the container's first element.
     * @param last I/O iterator to the container's end.
     * @param indices_first Input iterator to the first position.
     * @return
     */
    template<typename InputIt, typename OrderIt>
    constexpr inline auto order_from(InputIt first, InputIt last, OrderIt indices_first) {
        using std::next;
        const auto N = std::distance(first, last);
        const std::vector<typename std::iterator_traits<InputIt>::value_type> elems(first, last);
        for (auto i = 0; i < N; i++) { *first++ = elems[size_t(*next(indices_first, i))]; }
    }

    /**
     * Reorders a container taking elements from the given positions.
     * @param first I/O iterator to the container's first element.
     * @param N Container's size.
     * @param indices_first Input iterator to the first position.
     * @return
     */
    template<typename InputIt, typename OrderIt>
    constexpr inline auto order_from_n(InputIt first, size_t N, OrderIt indices_first) {
        return order_from(first, snext(first, N), indices_first);
    }

    // TODO: there must be a better way
    /**
     * Computes the ranks of a container's element.
     * @param first Input iterator to the container's first element.
     * @param last Input iterator to the container's end.
     * @param rank_first Output iterator to the ranks container.
     * @param compare Criterium to rank the elements. Must be usable by std::sort.
     */
    template<typename InputIt, typename RankIt, typename Compare>
    constexpr void rank(InputIt first, InputIt last, RankIt rank_first, Compare compare) {
        using std::distance, std::next;
        const auto N = distance(first, last);
        using index_type = typename std::iterator_traits<RankIt>::value_type;
        std::vector<index_type> argsort_indices((size_t(N)));
        argsort(first, last, argsort_indices.begin(), compare);
        std::iota(rank_first, next(rank_first, N), 0);
        order_to(rank_first, next(rank_first, N), argsort_indices.cbegin());
    }

    /**
     * Computes the ranks of a container's element.
     * @param first Input iterator to the container's first element.
     * @param N Container's size.
     * @param rank_first Output iterator to the ranks container.
     * @param compare Criterium to rank the elements. Must be usable by std::sort.
     */
    template<typename InputIt, typename RankIt, typename Compare>
    constexpr inline auto rank_n(InputIt first, size_t N, RankIt rank_first, Compare compare) {
        return rank(first, snext(first, N), rank_first, compare);
    }

    /**
     * Computes the ranks of a container's element sorting them with std::less.
     * @param first Input iterator to the container's first element.
     * @param last Input iterator to the container's end.
     * @param rank_first Output iterator to the ranks container.
     */
    template<typename InputIt, typename RankIt>
    constexpr inline auto rank(InputIt first, InputIt last, RankIt rank_first) {
        return rank(first, last, rank_first, std::less<>());
    }

    /**
     * Computes the ranks of a container's element sorting them with std::less.
     * @param first Input iterator to the container's first element.
     * @param N Container's size.
     * @param rank_first Output iterator to the ranks container.
     */
    template<typename InputIt, typename RankIt>
    constexpr inline auto rank_n(InputIt first, size_t N, RankIt rank_first) {
        return rank(first, snext(first, N), rank_first);
    }

    /**
     * Sorts the elements of
     * @tparam InputIt1
     * @tparam InputIt2
     * @tparam Compare
     * @param first_1
     * @param N
     * @param first_2
     * @param compare
     */
    template<typename InputIt1, typename InputIt2, typename Compare>
    void swap_order_by_rank_n(InputIt1 first_1, size_t N, InputIt2 first_2, Compare compare) {
        using std::next, std::distance;
        // Storing the original rank and sorting by that
        std::vector<int> rank_1(N);
        rank_n(first_1, N, rank_1.begin(), compare);
        order_to_n(first_1, N, rank_1.cbegin());
        std::vector<int> rank_2(N);
        rank_n(first_2, N, rank_2.begin(), compare);
        order_to_n(first_2, N, rank_2.cbegin());
        // Placing the elements in the order the others were:
        // the highest ranking element of container_1 is placed were the highest ranking element of
        // container_2 was, than the second-highest etc.
        // The sorting performed above was done in order to have the highest ranking element in first
        // position, the second-highest in second etc.
        order_from_n(first_1, N, rank_2.cbegin());
        order_from_n(first_2, N, rank_1.cbegin());
    }

    // TODO: there must be a better way
    template<typename InputIt1, typename InputIt2, typename Compare>
    void swap_order_by_rank(InputIt1 first_1, InputIt1 last_1, InputIt2 first_2, Compare compare) {
        return swap_order_by_rank_n(first_1, size_t(std::distance(first_1, last_1)), first_2,
                                    compare);
    }

    template<typename InputIt1, typename InputIt2>
    inline void swap_order_by_rank_n(InputIt1 first_1, size_t N, InputIt2 first_2) {
        return swap_order_by_rank_n(first_1, N, first_2, std::less<>());
    }

    template<typename InputIt1, typename InputIt2>
    inline void swap_order_by_rank(InputIt1 first_1, InputIt1 last_1, InputIt2 first_2) {
        swap_order_by_rank(first_1, last_1, first_2, std::less<>());
    }

    template<typename Val, typename It>
    inline bool in(const Val &x, It first, It last) {
        return first != last && std::find(first, last, x) != last;
    }


    namespace detail {
        template<typename ParentIt, typename ChildIt>
        constexpr inline void order_slice_as(ParentIt begin_parent, ParentIt end_parent,
                                             ParentIt begin_slice, ParentIt end_slice,
                                             ParentIt begin_index, ChildIt begin_child) {
            const auto i_size = std::distance(begin_parent, end_parent);
            auto copy_head = std::copy(begin_parent, begin_slice, begin_child);
            copy_head = std::copy_if(
                    begin_index, std::next(begin_index, i_size), copy_head,
                    [&](const auto &elem) { return in(elem, begin_slice, end_slice); });
            std::copy(end_slice, end_parent, copy_head);
        }
    }// namespace detail

    template<typename ParentIt, typename ChildIt>
    constexpr inline void cut_and_mix(ParentIt begin_parent1, ParentIt begin_parent2,
                                      ChildIt begin_child1, ChildIt begin_child2,
                                      const std::ptrdiff_t individual_size,
                                      std::ptrdiff_t start_slice1, std::ptrdiff_t start_slice2,
                                      std::ptrdiff_t slices_length) {
        using std::next;
        const auto begin_slice1 = next(begin_parent1, start_slice1);
        const auto end_slice1 = next(begin_slice1, slices_length);
        const auto end_parent1 = next(begin_parent1, individual_size);
        const auto begin_slice2 = next(begin_parent2, start_slice2);
        const auto end_slice2 = next(begin_slice2, slices_length);
        const auto end_parent2 = next(begin_parent2, individual_size);

        detail::order_slice_as(begin_parent1, end_parent1, begin_slice1, end_slice1, begin_parent2,
                               begin_child1);
        detail::order_slice_as(begin_parent2, end_parent2, begin_slice2, end_slice2, begin_parent1,
                               begin_child2);
    }

    namespace detail {
        template<typename ParentIt, typename ChildIt>
        inline ChildIt copy_slice_to(ParentIt begin_parent,
                                     const std::array<std::ptrdiff_t, 2> &slice_edges,
                                     ChildIt begin_child) {
            using std::next;
            return std::copy(next(begin_parent, slice_edges[0]), next(begin_parent, slice_edges[1]),
                             next(begin_child, slice_edges[0]));
        }

        template<typename ParentIt, typename ChildIt, class Condition>
        inline auto copy_n_if(ParentIt begin_from, ParentIt end_from, ChildIt begin_to, size_t n,
                              Condition condition) {
            while (n >= 1 && begin_from != end_from) {
                if (condition(*begin_from)) {
                    *begin_to++ = *begin_from;
                    n--;
                }
                begin_from++;
            }
            return std::make_pair(begin_from, begin_to);
        }

        template<typename Parent1It, typename Parent2It, typename ChildIt>
        inline void copy_import_slice(Parent1It begin_parent1, Parent1It end_parent1,
                                      Parent2It begin_slice_from,
                                      const std::ptrdiff_t start_slice_to,
                                      const std::ptrdiff_t slices_length, ChildIt begin_child) {
            const auto is_not_in_slice = [&](const auto &x) {
                return not in(x, begin_slice_from, std::next(begin_slice_from, slices_length));
            };
            auto [from_iter, to_iter] = detail::copy_n_if(begin_parent1, end_parent1, begin_child,
                                                          size_t(start_slice_to), is_not_in_slice);
            to_iter = std::copy_n(begin_slice_from, slices_length, to_iter);
            std::copy_if(from_iter, end_parent1, to_iter, is_not_in_slice);
        }
    }// namespace detail


    template<typename ParentIt, typename ChildIt>
    void exchange_slices(ParentIt begin_parent1, ParentIt begin_parent2, ChildIt begin_child1,
                         ChildIt begin_child2, const std::ptrdiff_t i_size,
                         std::ptrdiff_t start_slice1, std::ptrdiff_t start_slice2,
                         std::ptrdiff_t slices_length) {
        using std::next;
        const auto begin_slice1 = next(begin_parent1, start_slice1);
        const auto end_slice1 = next(begin_slice1, slices_length);
        const auto end_parent1 = next(begin_parent1, i_size);
        const auto begin_slice2 = next(begin_parent2, start_slice2);
        const auto end_slice2 = next(begin_slice2, slices_length);
        const auto end_parent2 = next(begin_parent2, i_size);

        detail::copy_import_slice(begin_parent1, end_parent1, begin_slice2, start_slice1,
                                  slices_length, begin_child1);
        detail::copy_import_slice(begin_parent2, end_parent2, begin_slice1, start_slice2,
                                  slices_length, begin_child2);
    }
}// namespace utils
#endif// GENETIC_TSP_UTILS_HPP
