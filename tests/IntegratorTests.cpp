//
// Created by Davide Nicoli on 21/10/22.
//
#include <limits>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "ariel_random/ariel_random.hpp"
#include "config.hpp"
#include "mc_integrators/integrator.hpp"
#include "samplers/MCMC/metropolis.hpp"
#include "transitions/uniform.hpp"

using namespace samplers::mcmc;
using namespace transitions;
using namespace mc_integrator;
using field = double;
using prob_space = double;

struct uniform_pdf {
    typedef field StateSpace;
    typedef field prob_space;
    explicit uniform_pdf(field limit = 1.0) : m_limit(limit){};
    [[nodiscard]] field logp(field x) const {
        if (x < 0.0 || x >= m_limit) return -std::numeric_limits<field>::infinity();
        return m_mlogvol;
    }
    field m_limit;
    field m_mlogvol{-std::log(m_limit)};
};

TEST_CASE("Integrator") {
    std::random_device rd;
    std::uniform_int_distribution<size_t> seed(0, 32000);
    SECTION("Metropolis sampling") {
        SECTION("Uniform") {
            ARandom rng(SEEDS_PATH "seed.in", PRIMES_PATH "primes32001.in", seed(rd));
            Metropolis sampler(field(0.5), uniform_pdf(),
                               UniformNear<prob_space, field>(/*radius=*/1.0));
            sampler.warmup(10000, rng);
            Integrator<field, decltype(sampler)> I(std::move(sampler));
            const auto [result, uncert] =
                    I([](const field & /*x*/) { return static_cast<field>(1.0); }, 1000, 10, rng);
            REQUIRE(result == Catch::Approx(1.0));
        }
        SECTION("Linear") {
            ARandom rng(SEEDS_PATH "seed.in", PRIMES_PATH "primes32001.in", seed(rd));
            Metropolis sampler(field(0.5), uniform_pdf(),
                               UniformNear<prob_space, field>(/*radius=*/0.5));
            Integrator<field, decltype(sampler)> I(std::move(sampler));
            const auto [result, n_blocks] =
                    I.integrate_to([](const field &x) { return x; }, 0.001, 100, rng);
            REQUIRE(std::get<0>(result) == Catch::Approx(0.5).epsilon(0.01));
        }
        SECTION("Sin") {
            ARandom rng(SEEDS_PATH "seed.in", PRIMES_PATH "primes32001.in", seed(rd));
            Metropolis sampler(field(M_PI / 2), uniform_pdf(M_PI),
                               UniformNear<prob_space, field>(/*radius=*/0.5));
            Integrator<field, decltype(sampler)> I(std::move(sampler));
            const auto [result, n_blocks] = I.integrate_to(
                    [](const field &x) { return std::sin(x) * M_PI; }, 0.001, 10, rng);
            REQUIRE(std::get<0>(result) == Catch::Approx(2.0).epsilon(0.01));
        }
        SECTION("Cos") {
            ARandom rng(SEEDS_PATH "seed.in", PRIMES_PATH "primes32001.in", seed(rd));
            Metropolis sampler(field(M_PI / 2), uniform_pdf(M_PI),
                               UniformNear<prob_space, field>(/*radius=*/1.0));
            Integrator<field, decltype(sampler)> I(std::move(sampler));
            const auto [result, n_blocks] = I.integrate_to(
                    [](const field &x) { return std::cos(x) * M_PI; }, 0.005, 10, rng);
            // Catch::Approx has a weird behavior when applied to 0
            CHECK(std::get<0>(result) + 1 == Catch::Approx(1.0).epsilon(0.1));
        }
    }
}