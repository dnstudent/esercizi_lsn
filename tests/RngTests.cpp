//
// Created by Davide Nicoli on 27/07/22.
//
#include <fstream>
#include <random>

#include <catch2/catch_test_macros.hpp>

#include "ariel_random/ariel_random.hpp"
#include "config.hpp"
#include "distributions/uniform_int.hpp"

TEST_CASE("ARandom", "[rng]") {
    SECTION("normal") {
        ARandom rng(SEEDS_PATH "seed.in", PRIMES_PATH "Primes", 1);
        std::normal_distribution<double> gauss(0, 1);
        std::ofstream results("gauss.out");
        for (size_t i = 0; i < 10'000; i++) { results << gauss(rng) << std::endl; }
        results.close();
    }
    SECTION("uniform int") {
        ARandom rng(SEEDS_PATH "seed.in", PRIMES_PATH "Primes", 3);
        std::uniform_int_distribution<size_t> unif(0, 100);
        std::ofstream results("unif_int.out");
        for (size_t i = 0; i < 10'000; i++) { results << unif(rng) << '\n'; }
        results.close();
    }
    SECTION("uniform double") {
        ARandom rng(SEEDS_PATH "seed.in", PRIMES_PATH "Primes", 0);
        std::uniform_real_distribution<double> unif(0, 100);
        std::ofstream results("unif_double.out");
        for (size_t i = 0; i < 10'000; i++) { results << unif(rng) << '\n'; }
        results.close();
    }
    SECTION("bernoulli") {
        ARandom rng(SEEDS_PATH "seed.in", PRIMES_PATH "Primes", 0);
        std::bernoulli_distribution coin{};
        std::ofstream results("bernoulli.out");
        for (size_t i = 0; i < 10'000; i++) { results << coin(rng) << '\n'; }
        results.close();
    }
    SECTION("Read primes") {
        size_t a, b;
        read_primes(PRIMES_PATH "Primes", 0, a, b);
        CHECK(a == 2892);
        CHECK(b == 2587);
        read_primes(PRIMES_PATH "Primes", 3, a, b);
        CHECK(a == 2892);
        CHECK(b == 2899);
        read_primes(PRIMES_PATH "primes32001.in", 0, a, b);
        CHECK(a == 2896);
        CHECK(b == 1263);
        read_primes(PRIMES_PATH "primes32001.in", 4, a, b);
        CHECK(a == 2896);
        CHECK(b == 1221);
    }
    SECTION("Read seeds") {
        size_t a, b, c, d;
        read_seeds(SEEDS_PATH "seed.in", a, b, c, d);
        CHECK(a == 0);
        CHECK(b == 0);
        CHECK(c == 0);
        CHECK(d == 1);
    }
}