//
// Created by Davide Nicoli on 28/09/22.
//

#include <catch2/catch_test_macros.hpp>

#include <limits>
#include <type_traits>
#include <vector>


enum Variable { En, P, T, K };


template<typename field, Variable var>
struct var_out {
    typedef field type;
};

template<typename field>
struct var_out<field, K> {
    typedef std::vector<field> type;
};

template<typename field, Variable var>
using var_out_t = typename var_out<field, var>::type;

template<Variable... vars>
struct varlist {
    using type = varlist;
    static constexpr std::size_t size() noexcept { return sizeof...(vars); }
};

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

template<Variable var, typename VarList>
struct indexer {};

template<Variable var, Variable... oth>
struct indexer<var, varlist<oth...>> {
    static constexpr size_t value = find_idx<0, var, oth...>();
};

template<Variable var, typename VarList>
inline constexpr size_t indexer_v = indexer<var, VarList>::value;

TEST_CASE("Meta", "[meta]") {
    SECTION("indexer") {
        REQUIRE(find_idx<0, Variable::T, Variable::T>() == 0);
        REQUIRE(find_idx<0, Variable::T, Variable::P, Variable::P, Variable::T, Variable::En>() ==
                2);
        REQUIRE(find_idx<0, Variable::T, Variable::P>() == static_cast<size_t>(-1));
        REQUIRE(indexer_v<Variable::T, varlist<P, En, T>> == 2);
    }
}