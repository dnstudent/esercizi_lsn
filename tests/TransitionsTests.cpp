//
// Created by Davide Nicoli on 25/08/22.
//

#include <limits>
#include <random>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "ariel_random/ariel_random.hpp"
#include "config.hpp"
#include "transitions/gauss.hpp"
#include "transitions/uniform.hpp"

using namespace transitions;
typedef std::valarray<double> StateSpace;

TEST_CASE("Transitions", "[trans]") {
    //    std::random_device rnd;
    //    std::uniform_int_distribution<uint_fast32_t> seed(0, 350);
    ARandom rng1(SEEDS_PATH "seed.in", PRIMES_PATH "Primes", 1);
    std::minstd_rand rng2(1);// NOLINT(cert-msc51-cpp)
    std::valarray<double> a{.1, .1};
    std::valarray<double> b{.2, .3};

    SECTION("Uniform") {
        SECTION("logp") {
            SECTION("out") {
                const double r = .05;
                UniformNear<double, StateSpace> t(r, a.size());
                auto logp = t.logp(b, a);
                REQUIRE((std::isinf(logp) && logp < 0));
            }
            SECTION("part_out") {
                const double r = .11;
                UniformNear<double, StateSpace> t(r, a.size());
                auto logp = t.logp(b, a);
                REQUIRE((std::isinf(logp) && logp < 0));
            }
            SECTION("in") {
                const double r = 0.5;
                UniformNear<double, StateSpace> t(r, a.size());
                auto logp = t.logp(b, a);
                REQUIRE(logp == Catch::Approx(0.0));
            }
        }
        SECTION("sample") {
            SECTION("ariel") {
                bool out = false;
                UniformNear<double, StateSpace> t(1.0, a.size());
                for (size_t i = 0; i < 1'000'000; i++) {
                    auto x = t.sample(a, rng1);
                    if (std::isinf(t.logp(x, a))) {
                        out = true;
                        break;
                    }
                }
                REQUIRE_FALSE(out);
            }
            SECTION("std") {
                bool out = false;
                UniformNear<double, StateSpace> t(1.0, a.size());
                for (size_t i = 0; i < 100'000; i++) {
                    auto x = t.sample(a, rng2);
                    if (std::isinf(t.logp(x, a))) {
                        out = true;
                        break;
                    }
                }
                REQUIRE_FALSE(out);
            }
        }
    }

    SECTION("Gauss") {
        SECTION("logp") {
            GaussNear<double, StateSpace> t1(1.0, a.size());
            REQUIRE(t1.logp(b, a) == Catch::Approx(-1.8628770664093453));

            GaussNear<double, StateSpace> t2(2.0, a.size());
            REQUIRE(t2.logp(b, a) == Catch::Approx(-3.2304214275292362));

            GaussNear<double, StateSpace> t3(3.0, a.size());
            REQUIRE(t3.logp(b, a) == Catch::Approx(-4.037879421523343));
        }
    }
}