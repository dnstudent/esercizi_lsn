//
// Created by Davide Nicoli on 22/07/22.
//

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "molecular_systems/data_types/vectors.hpp"

TEST_CASE("Vectors", "[structs]") {
    SECTION("-=") {
        Vectors<double> v1({1, 2}, {3, 4}, {5, 6});
        Vectors<double> v2({1, 2}, {4, 3}, {1, 1});
        v1 -= v2;
        REQUIRE(v1 == Vectors<double>({0, 0}, {-1, 1}, {4, 5}));
    }
    SECTION("*=") {
        Vectors<double> v1({1, 2}, {3, -4}, {-5, 6});
        double x = 3.4;
        v1 *= x;
        REQUIRE(v1 ==
                Vectors<double>({1 * 3.4, 2 * 3.4}, {3 * 3.4, -4 * 3.4}, {-5 * 3.4, 6 * 3.4}));
    }
    SECTION("=") {
        Vectors<double> v1({1, 2}, {3, 4}, {5, 6});
        Vectors<double> v2({1, -2}, {4, 3}, {-1, 1});
        REQUIRE(!(v1 == v2));
        v1 = v2;
        REQUIRE(v1 == v2);
    }
    SECTION("full_norm2") {
        Vectors<double> v1({1, 2}, {3, 4}, {5, 6});
        REQUIRE(v1.full_norm2() == Catch::Approx(91));
    }
}