#pragma once

#include <algorithm>
#include <array>
#include <filesystem>
#include <functional>
#include <iostream>
#include <iterator>
#include <numeric>
#include <type_traits>
#include <utility>
#include <vector>

#include <rapidcsv.h>

namespace csv = rapidcsv;
namespace fs = std::filesystem;

namespace utils {
    /**
     * Utility to do std::next of an unsigned size.
     */
    template<typename It>
    constexpr inline auto snext(It it, size_t N) {
        using diff_t = typename std::iterator_traits<It>::difference_type;
        return std::next(it, static_cast<diff_t>(N));
    }

    /**
     * Computes the sample average
     * @tparam var_space Numeric field of the result.
     * @param first Sample beginning.
     * @param last Sample's past-the-end iterator.
     * @return Sample average.
     */
    template<typename var_space, typename It>
    constexpr inline auto average(It first, It last) {
        return static_cast<var_space>(std::reduce(first, last)) /
               static_cast<var_space>(std::distance(first, last));
    }

    template<typename var_space, typename WeightIt, typename ElemIt>
    constexpr inline auto weighted_average(WeightIt first_weight, WeightIt last_weight,
                                           ElemIt first_elem) {
        return std::transform_reduce(first_weight, last_weight, first_elem, var_space(0)) /
               std::reduce(first_weight, last_weight);
    }

    /**
     * Computes the average of the squared rejects.
     * @tparam out_space Type of the output.
     * @param first Sample beginning.
     * @param last Sample's past-the-end iterator
     * @return Average of the squared rejects.
     */
    template<typename out_space, typename It>
    constexpr auto variance(It first, It last) {
        using in_space = typename std::iterator_traits<It>::value_type;
        const size_t N = unsigned(std::distance(first, last));
        std::vector<in_space> x2(N);
        std::transform(first, last, x2.begin(), [](const auto xi) { return xi * xi; });
        return average<out_space>(x2.cbegin(), x2.cend()) -
               std::pow(average<out_space>(first, last), 2);
    }

    /**
     * Computes the variance of a sample given its average.
     * @tparam out_space Type of the output.
     * @param first Sample beginning.
     * @param last Sample's past-the-end iterator
     * @return Average of the squared rejects.
     */
    template<typename out_space, typename It>
    constexpr auto variance(It first, It last, out_space sample_avg) {
        using in_space = typename std::iterator_traits<It>::value_type;
        const size_t N = unsigned(std::distance(first, last));
        const auto s2_avg =
                out_space(std::transform_reduce(first, last, in_space(0), std::plus<>(),
                                                [](const auto xi) { return xi * xi; })) /
                out_space(N);
        return s2_avg - std::pow(sample_avg, 2);
    }

    /*
    template<typename InputIt, typename OutputIt, typename Compare>
    constexpr inline void argsort(InputIt first, InputIt last, OutputIt d_first, Compare comp) {
        const auto N = std::distance(first, last);
        const auto d_last = next(d_first, N);
        std::iota(d_first, d_last, 0);
        std::stable_sort(d_first, d_last, [&](const auto left, const auto right) -> bool {
            // sort indices according to corresponding array element
            return comp(*next(first, left), *next(first, right));
        });
    }

    template<typename InputIt, typename OutputIt>
    constexpr inline void argsort(InputIt first, InputIt last, OutputIt d_first) {
        argsort(first, last, d_first, std::less<>());
    }
    */
    /**
     * L2 norm squared
     * @param v Input vector.
     */
    template<typename Vector>
    constexpr inline typename Vector::value_type norm2(const Vector &v) {
        return std::transform_reduce(std::begin(v), std::end(v), std::begin(v),
                                     typename Vector::value_type(0));
    }

    /**
     * Euclidean distance squared.
     * @param from
     * @param to
     */
    template<typename Point>
    constexpr inline auto distance2(const Point &from, const Point &to) {
        return std::transform_reduce(std::begin(from), std::end(from), std::begin(to),
                                     typename Point::value_type(0), std::plus<>(),
                                     [](const auto x, const auto y) {
                                         const auto d = y - x;
                                         return d * d;
                                     });
    }


    namespace detail {
        template<size_t I, class Fn, class... TupleLike>
        constexpr inline void _tuple_apply(Fn f, TupleLike &...tuples) {

            if constexpr (((I < std::tuple_size_v<TupleLike>) &&...)) {
                f(std::get<I>(tuples)...);
                _tuple_apply<I + 1>(f, tuples...);
            }
        }

        template<size_t I, class Fn, class OutTupleLike, class... InTupleLike>
        constexpr inline void _tuple_transform(Fn f, OutTupleLike &outs, InTupleLike &...tuples) {
            if constexpr (I < std::tuple_size_v<OutTupleLike> &&
                          ((I < std::tuple_size_v<InTupleLike>) &&...)) {
                std::get<I>(outs) = f(std::get<I>(tuples)...);
                _tuple_transform<I + 1>(f, outs, tuples...);
            }
        }
    }// namespace detail

    /**
     * Applies the same operation to groups of elements having the same index in different tuples.
     * @param f Operation to apply.
     * @param tuples Group of tuples to apply the operation to. The number of tuples must be equal to the number of parameters accepted by f.
     */
    template<class Fn, class... TupleLike>
    constexpr inline void tuple_apply(Fn f, TupleLike &...tuples) {
        detail::_tuple_apply<0>(f, tuples...);
    }


    /**
     * Applies the same operation to groups of elements having the same index in different tuples and stores the results.
     * @param f Operation to apply.
     * @param outs Tuple-like structure in which results will be saved.
     * @param ins Group of tuples to apply the operation to. The number of tuples must be equal to the number of parameters accepted by f.
     */
    template<class Fn, class OutTupleLike, class... InTupleLike>
    constexpr inline void tuple_transform(Fn f, OutTupleLike &outs, InTupleLike &...ins) {
        detail::_tuple_transform<0>(f, outs, ins...);
    }

    /**
     * Flattens a rank 2 tuple
     * @tparam Elems Tuples.
     * @param tuple A tuple of tuples.
     * @return A rank 1 tuple.
     */
    template<class... Elems>
    constexpr inline auto tuple_flatten(const std::tuple<Elems...> &tuple) {
        return std::apply([&](Elems const &...e) { return std::tuple_cat(e...); }, tuple);
    }


    /**
     * Appends the given tuple of vectors to a rapidcsv Document.
     * @param table Rapidcsv table.
     * @param names Column names.
     * @param columns Tuple-like of vectors.
     */
    template<size_t I = 0, class Tuple>
    constexpr inline void
    AppendColumns(csv::Document &table,
                  const std::array<std::string, std::tuple_size_v<Tuple>> &names,
                  const Tuple &columns) {
        tuple_apply(
                [&](const auto &column, const auto &name) {
                    const size_t c = table.GetColumnCount();
                    table.SetColumn(c, column);
                    table.SetColumnName(c, name);
                },
                columns, names);
    }

    /**
     * Appends the given stats to a rapidcsv Document.
     * @param table Rapidcsv document.
     * @param var_names Names of the variables.
     * @param stat_names Names of the statistics.
     * @param measures_stats Tuple of tuples; the first axis changes the variable, the second changes the statistic.
     */
    template<size_t I = 0, size_t N_STATS, class MeasuresTuple>
    constexpr inline void
    AppendColumns(csv::Document &table,
                  const std::array<std::string, std::tuple_size_v<MeasuresTuple>> &var_names,
                  const std::array<std::string, N_STATS> &stat_names,
                  const MeasuresTuple &measures_stats) {
        tuple_apply(
                [&](const auto &var_stats, const auto &var_name) {
                    std::array<std::string, N_STATS> col_names;
                    for (size_t i = 0; i < N_STATS; i++) {
                        col_names[i] = var_name + "_" + stat_names[i];
                    }
                    AppendColumns(table, col_names, var_stats);
                },
                measures_stats, var_names);
    }


    /**
     * Push back each element of values to the corresponding vector in to.
     * @param values Values to append.
     * @param to Vectors where values will be appended.
     */
    template<size_t I = 0, typename SingleTuple, typename VectorTuple>
    constexpr inline void tuple_push_back(const SingleTuple &values, VectorTuple &to) {
        tuple_apply([](const auto &value, auto &v) { v.push_back(value); }, values, to);
    }

    /**
     * Requires the existence of path, else throws
     * @param path Path.
     */
    inline void require_existence(const fs::path &path) {
        if (!fs::exists(path)) {
            throw std::runtime_error("The directory: " + path.string() + " do not exist.");
        }
    }

    /**
     * Computes autocorrelation for a given lag using the provided sum of the rejects squared. It assumes
     * the sample contains the rejects.
     * @param first Sample beginning.
     * @param last Sample's past-the-end iterator.
     * @param sample_rej2_sum Sum of the rejects squared.
     * @param lag Lag for which autocorrelation is computed.
     * @return Autocorrelation for the given lag.
     */
    template<typename InputIt>
    constexpr inline auto
    _autocorrelation(InputIt first, InputIt last,
                     const typename std::iterator_traits<InputIt>::value_type sample_rej2_sum,
                     const size_t lag) {
        const auto dlag = typename std::iterator_traits<InputIt>::difference_type(lag);
        return std::transform_reduce(first, std::prev(last, dlag), std::next(first, dlag),
                                     typename std::iterator_traits<InputIt>::value_type(0)) /
               sample_rej2_sum;
    }

    /**
     * Computes the autocorrelation function of a given sample for the specified number of lags using
     * the definition. FFT would be way better.
     * @param first Sample beginning.
     * @param last Sample's past-the-end iterator.
     * @param first_out Output beginning.
     * @param n_lags Number of lags.
     */
    template<typename InputIt, typename OutputIt>
    constexpr inline void autocorrelation_fn(InputIt first, InputIt last, OutputIt first_out,
                                             size_t n_lags = 0) {
        if (n_lags == 0) n_lags = size_t(std::distance(first, last));
        using ovar_space = typename std::iterator_traits<OutputIt>::value_type;
        const auto t_max = size_t(std::distance(first, last));
        const auto sample_avg = average<ovar_space>(first, last);
        std::vector<ovar_space> y(t_max);
        std::transform(first, last, y.begin(),
                       [=](const auto xt) { return ovar_space(xt) - sample_avg; });
        const auto sample_rej2_sum =
                std::transform_reduce(y.cbegin(), y.cend(), y.cbegin(), ovar_space(0));
        for (size_t t = 0; t < n_lags; t++)
            *first_out++ = _autocorrelation(y.cbegin(), y.cend(), sample_rej2_sum, t);
    }

    /**
     * Computes the autocorrelation for every column of a given csv file and stores it in another csv file.
     * @param input_path Path to the series.
     * @param output_path Path to the output.
     * @param n_lags Number of lags.
     * @param skip Number of rows to skip from the beginning of input.
     */
    template<typename real>
    inline void autocorrelation_from(const fs::path &input_path, const fs::path &output_path,
                                     size_t n_lags, size_t skip) {
        csv::Document in(input_path);
        csv::Document out;
        for (size_t i = 0; i < in.GetColumnCount(); i++) {
            const auto data = in.template GetColumn<real>(i);
            std::vector<real> ac_fn(n_lags);
            autocorrelation_fn(snext(data.cbegin(), skip), data.cend(), ac_fn.begin(), n_lags);
            out.InsertColumn(i, std::move(ac_fn), in.GetColumnName(signed(i)));
        }
        out.RemoveColumn(out.GetColumnCount() - 1);
        if (!fs::exists(output_path.parent_path()))
            fs::create_directories(output_path.parent_path());
        out.Save(output_path);
    }


    /**
     * Every element in sample must be >= a and < b
     * @tparam XIt
     * @tparam BinsIt
     * @param sample_first
     * @param sample_last
     * @param bins_first
     * @param bins_last
     * @param a
     * @param b
     */
    template<typename XIt, typename BinsIt, typename EdgeIt>
    inline void histogram(XIt sample_first, XIt sample_last, BinsIt bins_first, BinsIt bins_last,
                          EdgeIt edge_first, typename XIt::value_type a,
                          typename XIt::value_type b) {
        using field = typename XIt::value_type;
        static_assert(std::is_floating_point_v<field>);
        using diff_t = typename BinsIt::difference_type;
        const auto n_bins = static_cast<size_t>(std::distance(bins_first, bins_last));
        std::fill(bins_first, bins_last, 0);
        const auto bin_size = (b - a) / static_cast<field>(n_bins);
        while (sample_first != sample_last) {
            *std::next(bins_first, static_cast<diff_t>((*sample_first - a) / bin_size)) += 1;
            sample_first++;
        }
        std::generate(edge_first, snext(edge_first, n_bins), [e = a, bin_size]() mutable {
            const auto e_prev = e;
            e += bin_size;
            return e_prev;
        });
    }


}// namespace utils
