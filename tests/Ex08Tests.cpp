//
// Created by Davide Nicoli on 22/10/22.
//

#include <algorithm>
#include <filesystem>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <rapidcsv.h>

#include "config.hpp"
#include "distributions/exercises.hpp"

using namespace distributions::ex08;
using namespace functions::ex08;
namespace csv = rapidcsv;
namespace fs = std::filesystem;

using field = double;
using prob_space = double;

TEST_CASE("Exercise 08") {
    SECTION("Integrand") {
        SECTION("1,1") {
            Integrand<field> integrand(1, 1);
            csv::Document table(fs::path(TESTS_PATH) / "ex08" / "integrand_1_1.csv");
            const auto xs = table.GetColumn<field>("x");
            std::vector<field> ys(xs.size());
            std::transform(xs.cbegin(), xs.cend(), ys.begin(), integrand);
            std::vector<field> expected = table.GetColumn<field>("y");
            for (size_t i = 0; i < expected.size(); i++)
                REQUIRE(ys[i] == Catch::Approx(expected[i]));
        }
        SECTION("2,1.5") {
            Integrand<field> integrand(2, 1.5);
            csv::Document table(fs::path(TESTS_PATH) / "ex08" / "integrand_2_1.5.csv");
            const auto xs = table.GetColumn<field>("x");
            std::vector<field> ys(xs.size());
            std::transform(xs.cbegin(), xs.cend(), ys.begin(), integrand);
            std::vector<field> expected = table.GetColumn<field>("y");
            for (size_t i = 0; i < expected.size(); i++)
                REQUIRE(ys[i] == Catch::Approx(expected[i]));
        }
    }
    SECTION("Pdf") {
        SECTION("1,1") {
            Trial<field> pdf(1, 1);
            csv::Document table(fs::path(TESTS_PATH) / "ex08" / "trial_1_1.csv");
            const auto xs = table.GetColumn<field>("x");
            std::vector<field> ys(xs.size());
            std::transform(xs.cbegin(), xs.cend(), ys.begin(),
                           [&](const auto x) { return pdf.logp(x); });
            std::vector<field> expected = table.GetColumn<field>("y");
            for (size_t i = 0; i < expected.size(); i++)
                CHECK(ys[i] == Catch::Approx(expected[i]).epsilon(0.00001));
        }
        SECTION("2,1.5") {
            Trial<field> pdf(2, 1.5);
            csv::Document table(fs::path(TESTS_PATH) / "ex08" / "trial_2_1.5.csv");
            const auto xs = table.GetColumn<field>("x");
            std::vector<field> ys(xs.size());
            std::transform(xs.cbegin(), xs.cend(), ys.begin(),
                           [&](const auto x) { return pdf.logp(x); });
            std::vector<field> expected = table.GetColumn<field>("y");
            for (size_t i = 0; i < expected.size(); i++)
                CHECK(ys[i] == Catch::Approx(expected[i]).epsilon(0.0001));
        }
        SECTION("2,0") {
            Trial<field> pdf(2, 0.0005);
            csv::Document table(fs::path(TESTS_PATH) / "ex08" / "trial_2_0.csv");
            const auto xs = table.GetColumn<field>("x");
            std::vector<field> ys(xs.size());
            std::transform(xs.cbegin(), xs.cend(), ys.begin(),
                           [&](const auto x) { return pdf.logp(x); });
            std::vector<field> expected = table.GetColumn<field>("y");
            for (size_t i = 0; i < expected.size(); i++)
                CHECK(ys[i] == Catch::Approx(expected[i]).epsilon(0.0001));
        }
    }
}