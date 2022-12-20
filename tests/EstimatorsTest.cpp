//
// Created by Davide Nicoli on 04/07/22.
//

#include <valarray>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "estimators/mean.hpp"


TEST_CASE("Testing estimators", "[estimators]") {
    const std::valarray<double> block1{.1, .2, .3, .4};
    const double av1 = block1.sum() / double(block1.size());
    const std::valarray<double> block2{.2, .3, .2, .3};
    const double av2 = block2.sum() / double(block2.size());
    const std::valarray<double> block3{.3, .3, .4, .5};
    const double av3 = block3.sum() / double(block3.size());
    const std::valarray<double> avs{av1, av2, av3};
    SECTION("Mean") {
        const auto [estimate, error] =
                estimators::Average<double>()(std::begin(block1), std::end(block1));
        REQUIRE(estimate == Catch::Approx(double(0.25)));
        REQUIRE(error == Catch::Approx(0.06454972243679029));
    }

    SECTION("BlockMean") {
        estimators::ProgAvg<double> est;
        auto [estimate, error] = est(std::begin(block1), std::end(block1));
        REQUIRE(estimate == Catch::Approx(double(0.25)));
        REQUIRE(error == Catch::Approx(0));
        std::tie(estimate, error) = est(std::begin(block2), std::end(block2));
        REQUIRE(estimate == Catch::Approx(double(0.25)));
        REQUIRE(error == Catch::Approx(0));
        std::tie(estimate, error) = est(std::begin(block3), std::end(block3));
        REQUIRE(estimate == Catch::Approx(double(0.2916666666666667)));
        REQUIRE(error == Catch::Approx(0.041666666666666595));

        const auto [mean_est, mean_err] =
                estimators::Average<double>()(std::begin(avs), std::end(avs));

        REQUIRE(mean_est == Catch::Approx(estimate));
        REQUIRE(mean_err == Catch::Approx(error));
    }
}