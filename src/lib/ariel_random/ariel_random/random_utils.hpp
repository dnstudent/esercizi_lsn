//
// Created by Davide Nicoli on 20/05/22.
//

#ifndef GENETIC_TSP_RANDOM_UTILS_HPP
#define GENETIC_TSP_RANDOM_UTILS_HPP

#include <cstddef>
#include <fstream>
#include <iterator>
#include <numeric>
#include <vector>

#include "ariel_random.hpp"

// Fisher-Yates algorithm with ARandom
template <typename It>
void random_shuffle(It first, const It last, ARandom &rng) {
  const auto N = std::distance(first, last);
  if (N <= 1)
    return;
  while (first != std::prev(last)) {
    const auto offset = typename std::iterator_traits<It>::difference_type(
        rng.Rannyu(0, double(std::distance(first, last))));
    std::iter_swap(first, std::next(first, offset));
    first++;
  }
}

template <typename uint>
[[maybe_unused]] uint read_primes(uint &a, uint &b,
                                  const std::string_view &primes_source,
                                  size_t primes_line = 0) {
  std::ifstream primes((std::string(primes_source)));
  if (primes.is_open()) {
    for (size_t i = 0; i + 1 < primes_line && !primes.eof(); i++)
      primes >> a >> b;
    if (!primes.eof()) {
      primes >> a >> b;
      primes.close();
      return b4096tob10(uint(0), uint(0), a, b);
    } else {
      throw std::runtime_error("The primes file was too short");
    }
  } else {
    try {
      primes.close();
    } catch (...) {
    }
    throw std::runtime_error("Could not open " + std::string(primes_source));
  }
}

template <typename uint>
[[maybe_unused]] uint read_seed(const std::string_view &seeds_source) {
  std::ifstream input((std::string(seeds_source)));
  std::string property;
  if (input.is_open()) {
    while (!input.eof()) {
      input >> property;
      if (property == "RANDOMSEED") {
        uint a, b, c, d;
        input >> a >> b >> c >> d;
        input.close();
        return b4096to10(a, b, c, d);
      }
    }
    throw std::runtime_error("Reched the end of " + std::string(seeds_source) +
                             " without finding the correct specification");
  } else {
    try {
      input.close();
    } catch (...) {
    }
    throw std::runtime_error("Could not open " + std::string(seeds_source));
  }
}

#include "ariel_random.hpp"

#endif // GENETIC_TSP_RANDOM_UTILS_HPP
