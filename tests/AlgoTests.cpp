//
// Created by Davide Nicoli on 22/07/22.
//

#include <cmath>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <rapidcsv.h>

#include "config.hpp"
#include "molecular_systems/algos.hpp"
#include "molecular_systems/data_types/settings.hpp"
#include "molecular_systems/data_types/vectors.hpp"
#include "utils.hpp"

TEST_CASE("Testing algorithms", "[algo]") {
    SECTION("pbc") {
        double box_edge{1};
        double x1{1.3};
        REQUIRE(PBC(x1, box_edge) == Catch::Approx(0.3));
        double x2{0.6};
        double x3{-0.5};
        REQUIRE(PBC(std::abs(x3 - x2), box_edge) == Catch::Approx(0.1));
    }
    SECTION("autocorrelation") {
        rapidcsv::Document flights(TESTS_PATH "flights.csv");
        const std::vector<size_t> pass = flights.GetColumn<size_t>("passengers");
        std::vector<double> autoc(pass.size());
        utils::autocorrelation_fn(pass.cbegin(), pass.cend(), autoc.begin());
        std::vector<double> expected{1.,         0.94804734, 0.87557484, 0.80668116,
                                     0.75262542, 0.71376997, 0.6817336,  0.66290439,
                                     0.65561048, 0.67094833, 0.70271992};
        for (size_t i = 0; i < expected.size(); i++)
            REQUIRE(autoc[i] == Catch::Approx(expected[i]));
    }
}

TEST_CASE("Physics", "[algo]") {
    SECTION("LJ") {
        SimulationSettings<double> ss(2, 1, 1, 5.0, 0.2, 0.05 /*, true*/);
        Vectors<double> positions1{{0, 0.1 * 12}, {0, 0}, {0, 0}};
        Vectors<double> positions2{{0, 0}, {0, 0.1 * 12}, {0, 0}};
        Vectors<double> positions3{{0, 0}, {0, 0}, {0, 0.1 * 12}};
        auto f1 = [&](const auto &p) { return LJ_potential<false>(0, {0, 0, 0}, p, ss); };
        auto f2 = [&](const auto &p) { return LJ_potential<true>(0, {0, 0, 0}, p, ss); };
        REQUIRE(f1(positions1) == Catch::Approx(f1(positions2)));
        REQUIRE(f1(positions2) == Catch::Approx(f1(positions3)));
        REQUIRE(f1(positions1) == Catch::Approx(f1(positions3)));
        REQUIRE(f1(positions1) == Catch::Approx(-0.8909652875830756));
        REQUIRE(f2(positions1) == Catch::Approx(f2(positions2)));
        REQUIRE(f2(positions2) == Catch::Approx(f2(positions3)));
        REQUIRE(f2(positions1) == Catch::Approx(f2(positions3)));
    }
}