//
// Created by Davide Nicoli on 24/07/22.
//

#include <catch2/catch_test_macros.hpp>

#include "molecular_systems/algos.hpp"
#include "molecular_systems/data_types/vectors.hpp"

TEST_CASE("MD", "[md]") {
    SECTION("utils") {
        SECTION("Previous positions") {
            Vectors<double> positions{{0.0, 0.1}, {0.0, 0.1}, {0.0, 0.1}};
            Vectors<double> velocities{{0.1, 0.1}, {0.1, 0.1}, {0.1, 0.1}};
            Vectors<double> pp(2);
            double dt = 0.1;
            double box_edge = 1;
            compute_previous_positions(positions, velocities, dt, box_edge, pp);
        }
    }
}
