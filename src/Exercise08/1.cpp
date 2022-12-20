//
// Created by Davide Nicoli on 19/10/22.
//
#include <iostream>
#include <random>

#include "ariel_random/ariel_random.hpp"
#include "config.hpp"
#include "distributions/exercises.hpp"
#include "mc_integrators/integrator.hpp"
#include "samplers/MCMC/metropolis.hpp"
#include "transitions/uniform.hpp"

using namespace samplers::mcmc;
using namespace distributions::ex08;
using namespace functions::ex08;
using namespace transitions;
using namespace mc_integrator;

typedef double field;
typedef double prob_space;

int main(int /*argc*/, char ** /*argv[]*/) {
    std::random_device rd;
    std::uniform_int_distribution<size_t> seed(0, 32000);
    ARandom rng(SEEDS_PATH "seed.in", PRIMES_PATH "primes32001.in", seed(rd));
    Integrand<field> Hpsi(/*mu=*/1.1, /*sigma=*/1);
    Metropolis sampler(0, Trial<field>(/*mu=*/1, /*sigma=*/1),
                       UniformNear<prob_space, field>(/*radius=*/1));
    sampler.warmup(10000, rng);
    Integrator<field, decltype(sampler)> I(std::move(sampler));
    const auto result = I(Hpsi, 10000, 1000, rng);
    std::cout << "Estimate: " << std::get<0>(result) << "\nUncertainty: " << std::get<1>(result)
              << std::endl;
    return 0;
}