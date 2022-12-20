//
// Created by Davide Nicoli on 09/09/22.
//

#include <random>

#include "config.hpp"
#include "estimators/mean.hpp"
#include "models/ising/1D/ising.hpp"
#include "models/ising/1D/simulation.hpp"
#include "models/ising/1D/variables.hpp"
#include "samplers/MCMC/metropolis.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

using namespace ising;
using namespace ising::D1;
using namespace spins;
using namespace samplers::mcmc;

using Ising1D = Ising<double, 1>;
using M_var = Magnetization<double, ProgAvg<double>>;

TEST_CASE("Ising", "[models]") {
    std::minstd_rand rng;
    Ising1D i1(3, down, 1, 0.1, 1);
    Ising1D i2({down, down, up}, 1, 0.1, 1);
    Ising1D i3({up, down, down}, 1, 0.1, 1);
    SECTION("∆E") {
        REQUIRE(i1.flip_dE(2) == Catch::Approx(i2.energy() - i1.energy()));
        REQUIRE(i1.flip_dE(0) == Catch::Approx(i3.energy() - i1.energy()));
    }
    SECTION("E") {
        Ising1D i4(4, down, 1, 0, 1);
        REQUIRE(i4.energy() == Catch::Approx(-4));
        Ising1D i5({down, down, down, up, up, up, down}, 1, 0, 1);
        REQUIRE(i5.energy() == Catch::Approx(-3));
        REQUIRE(i1.energy() == Catch::Approx(-2.7));
        REQUIRE(i2.energy() == Catch::Approx(i3.energy()));
    }
    SECTION("Spin sum") {
        REQUIRE(D1::spin_sum(i1) == -3);
        REQUIRE(D1::spin_sum(i2) == -1);
        REQUIRE(D1::spin_sum(i2) == D1::spin_sum(i3));
    }
    SECTION("M") {
        auto i6 = std::make_shared<Ising1D>(5, down, 1, 0, 0);
        REQUIRE(D1::spin_sum(*i6) == -5);
        Simulator<false, true, false, double, SystemMetropolis<double>, M_var> sim(10, i6,
                                                                                   M_var(*i6));
        sim.run(1, 0, rng);
        sim.save_results(fs::path(RESULTS_DIR "test.csv"));
        sim.save_state(fs::path(RESULTS_DIR "state.csv"));
    }
}