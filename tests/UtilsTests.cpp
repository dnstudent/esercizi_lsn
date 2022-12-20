//
// Created by Davide Nicoli on 06/10/22.
//

#include <tuple>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <rapidcsv.h>

#include "utils.hpp"

template<typename x>
using v = std::vector<x>;
namespace csv = rapidcsv;

TEST_CASE("Utils", "[utils]") {
    SECTION("Tuple apply") {
        const auto a = std::make_tuple<int, bool, char, std::string>(1, true, 'a', "str");
        auto av = std::make_tuple<v<int>, v<bool>, v<char>, v<std::string>>({}, {}, {}, {});
        utils::tuple_push_back(a, av);
        REQUIRE(av == decltype(av){{1}, {true}, {'a'}, {"str"}});
    }

    SECTION("Column append") {
        csv::Document table;
        auto cols = std::make_tuple(v<int>{1, 2, 3}, v<char>{'b', 'c', 'a'},
                                    v<float>{3.4f, 22.0f, 1.0f});
        utils::AppendColumns(table, {"ints", "chars", "floats"}, cols);
        REQUIRE(true == true);
    }

    SECTION("Histogram") {
        SECTION("Zero") {
            const std::vector<double> sample{};
            std::vector<double> edges(3);
            std::vector<size_t> bins(3, 0);
            utils::histogram(sample.begin(), sample.end(), bins.begin(), bins.end(), edges.begin(),
                             double(0), double(6));
            CHECK(bins == std::vector<size_t>{0, 0, 0});
        }
        SECTION("Simple") {
            const std::vector<double> sample{0, 1, 2, 3, 4, 5, 2, 3, 1, 5};
            std::vector<double> edges(3);
            std::vector<size_t> bins(3, 0);
            utils::histogram(sample.begin(), sample.end(), bins.begin(), bins.end(), edges.begin(),
                             double(0), double(6));
            CHECK(bins == std::vector<size_t>{3, 4, 3});
        }
    }
}