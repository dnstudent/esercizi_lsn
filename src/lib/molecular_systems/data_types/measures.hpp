//
// Created by Davide Nicoli on 28/09/22.
//

#ifndef ESERCIZI_LSN_MS_MEASURES_HPP
#define ESERCIZI_LSN_MS_MEASURES_HPP

#include <cstddef>
#include <vector>

#include <rapidcsv.h>

#include "utils.hpp"

namespace csv = rapidcsv;

enum Variable {
    PotentialEnergy,
    KineticEnergy,
    TotalEnergy,
    Pressure,
    Temperature,
    RadialFn,
};


template<typename field, Variable var>
struct var_out {
    typedef field type;
};

template<typename field>
struct var_out<field, RadialFn> {
    typedef std::vector<field> type;
};

template<typename field, Variable var>
using var_out_t = typename var_out<field, var>::type;

template<Variable... vars>
struct varlist {
    using type = varlist;
    static constexpr std::size_t size() noexcept { return sizeof...(vars); }
};

namespace detail {
    template<size_t N = 0, Variable value>
    static constexpr size_t find_idx() {
        return static_cast<size_t>(-1);
    }

    template<size_t N = 0, Variable value, Variable first, Variable... oth>
    static constexpr size_t find_idx() {
        if constexpr (value == first) return N;
        else
            return find_idx<N + 1, value, oth...>();
    }
}// namespace detail

template<Variable var, typename VarList>
struct indexer {};

template<Variable var, Variable... oth>
struct indexer<var, varlist<oth...>> {
    static constexpr size_t value = detail::find_idx<0, var, oth...>();
};

template<Variable var, typename VarList>
inline constexpr size_t indexer_v = indexer<var, VarList>::value;

template<typename field, Variable... vars>
class MeasureOutputs {
    typedef varlist<vars...> VarList;

public:
    MeasureOutputs() = default;
    explicit MeasureOutputs(size_t n_measures) {
        utils::tuple_apply([&](auto &ms) { ms.reserve(n_measures); }, m_measures);
    }

    constexpr void clear() {
        utils::tuple_apply([&](auto &measures) { measures.clear(); }, m_measures);
    }

    template<Variable var>
    constexpr void push_measure(const var_out_t<field, var> &measure) {
        if constexpr (has_member<var>()) {
            std::get<indexer_v<var, VarList>>(m_measures).push_back(measure);
        }
    }

    template<Variable var>
    constexpr void push_measure(var_out_t<field, var> &&measure) {
        if constexpr (has_member<var>()) {
            std::get<indexer_v<var, VarList>>(m_measures)
                    .push_back(std::forward<var_out_t<field, var>>(measure));
        }
    }

    // NOTE: this is terrible
    template<Variable var>
    constexpr void push_vector(const std::vector<field> &measure) {
        if constexpr (has_member<var>()) {
            for (size_t bin_idx = 0; bin_idx < measure.size(); bin_idx++) {
                std::get<indexer_v<var, VarList>>(m_measures)[bin_idx].push_back(measure[bin_idx]);
            }
        }
    }
    template<Variable var>
    constexpr void init_vector(size_t n_bins) {
        if constexpr (has_member<var>()) {
            std::get<indexer_v<var, VarList>>(m_measures).resize(n_bins);
        }
    }

    template<Variable var>
    static constexpr bool has_member() {
        return ((var == vars) || ...);
    }

    template<Variable var>
    constexpr auto &get_measures() const {
        static_assert(has_member<var>());
        return std::get<indexer_v<var, VarList>>(m_measures);
    }

    auto &all_measures() const { return m_measures; }

    static constexpr size_t N_VARS = sizeof...(vars);

    // Quick and dirty trick to deduce the number of estimators in other classes
    static constexpr size_t N_SCALARS =
            N_VARS - static_cast<size_t>(has_member<Variable::RadialFn>());

private:
    std::tuple<std::vector<var_out_t<field, vars>>...> m_measures;
};

#endif//ESERCIZI_LSN_MS_MEASURES_HPP
