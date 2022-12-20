//
// Created by Davide Nicoli on 10/06/22.
//

#include <algorithm>
#include <array>
#include <iostream>
#include <random>
#include <valarray>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "config.hpp"

#include "ariel_random/ariel_random.hpp"
#include "genetic/genetic_utils.hpp"
#include "genetic_tsp/crossovers.hpp"
#include "genetic_tsp/tsp_ga.hpp"

typedef std::valarray<float> point;
typedef std::array<point, 50> Coordinates;
using namespace utils;

#define GEN_IND(name, rng, Individual)                                                             \
    Individual name;                                                                               \
    std::iota((name).begin(), (name).end(), 1);                                                    \
    std::shuffle((name).begin(), (name).end(), (rng));

template<typename Individual>
bool is_valid_individual(Individual &i) {
    std::sort(i.begin(), i.end());
    Individual base{};
    std::iota(base.begin(), base.end(), 1);
    return i == base;
}

TEST_CASE("Genetic algorithm", "[ga]") {
    SECTION("Distance") {
        std::array<std::array<float, 2>, 5> coordinates{{{0, 0}, {0, 1}, {1, 1}, {1, 0}, {0, 0}}};
        TSP ga(coordinates, MyCrossover2<5>{});
        using Individual = decltype(ga)::Individual;
        Individual square_path{1, 2, 3, 4};
        REQUIRE(ga.fitness(square_path) == Catch::Approx(0.25));
        Individual crossed_path{1, 3, 2, 4};
        REQUIRE(ga.fitness(crossed_path) == Catch::Approx(1.0 / (1.0 + M_SQRT2 + 1.0 + M_SQRT2)));
    }
    SECTION("Crossover") {
        // ARandom rng(SEEDS_PATH "seed.in", PRIMES_PATH "primes32001.in", 5U);
        std::random_device rnd;
        std::minstd_rand rng(rnd());
        std::uniform_real_distribution<float> p_dist(0, 1);
        Coordinates coordinates{};
        std::generate(coordinates.begin(), coordinates.end(), [&]() {
            return point{p_dist(rng), p_dist(rng)};
        });
        TSP ga(coordinates, MyCrossover2<50>{});
        using Individual = decltype(ga)::Individual;
        for (auto i = 0; i < 1'000'000; i++) {
            GEN_IND(first, rng, Individual)
            GEN_IND(second, rng, Individual)
            Individual first_child;
            Individual second_child;
            REQUIRE_NOTHROW(ga.crossover(first, second, first_child, second_child, rng));
            REQUIRE(is_valid_individual(first_child));
            REQUIRE(is_valid_individual(second_child));
        }
    }
}
TEST_CASE("Genetic utils", "[gu]") {
    SECTION("Cut and mix") {
        std::vector<uint8_t> parent_1{1, 2, 3, 4, 5, 6, 7};
        std::vector<uint8_t> parent_2{2, 6, 3, 1, 4, 5, 7};
        std::vector<uint8_t> child_1(parent_1.size());
        std::vector<uint8_t> child_2(parent_2.size());

        cut_and_mix(parent_1.cbegin(), parent_2.cbegin(), child_1.begin(), child_2.begin(),
                    static_cast<ptrdiff_t>(parent_1.size()), 0, 0, 7);

        REQUIRE(child_1 == parent_2);
        REQUIRE(child_2 == parent_1);

        cut_and_mix(parent_1.cbegin(), parent_2.cbegin(), child_1.begin(), child_2.begin(),
                    static_cast<ptrdiff_t>(parent_1.size()), 4, 4, 3);

        REQUIRE(child_1 == std::vector<uint8_t>{1, 2, 3, 4, 6, 5, 7});
        REQUIRE(child_2 == std::vector<uint8_t>{2, 6, 3, 1, 4, 5, 7});
    }
    SECTION("copy_n_if") {
        std::vector<int> a{1, 3, 5, 2, 4, 6, 1, 2, 6, 1, 3, 1, 3, 1, 2, 3, 1, 2, 4, 5, 1, 2, 3, 1};
        std::vector<int> b(4);
        REQUIRE_NOTHROW(utils::detail::copy_n_if(a.begin(), a.end(), b.begin(), 4,
                                                 [](const auto x) { return x < 4; }));
        REQUIRE(b == std::vector<int>{1, 3, 2, 1});
    }
}