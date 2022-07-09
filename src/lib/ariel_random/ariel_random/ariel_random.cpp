/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include "ariel_random.hpp"

#define twom12 0.000244140625

/**
 * Initializes the number generator with the data in the provided files.
 * @param seeds_source Path to a file containing the seed. Must be on a line starting with "RANDOMSEED" followed by four base-2^12 digits separated by spaces.
 * @param primes_source Path to a file containing primes. Each line of the file must list the two base-2^12 digits representing a prime < 2^24, separated by a single space.
 * @param primes_line The line of the primes file to use, starting from 1. If 0 use the default prime (11848219)
 */
ARandom::ARandom(const std::string_view &seeds_source, const std::string_view &primes_source,
                 size_t primes_line) {
    std::ifstream primes((std::string(primes_source)));
    if (primes.is_open()) {
        result_type bucket;
        for (size_t i = 0; i < primes_line && !primes.eof(); i++) primes >> bucket >> bucket;
        if (!primes.eof()) {
            primes >> m_p3 >> m_p4;
            std::cout << "Prime " << b4096tob10<result_type>(0, 0, m_p3, m_p4) << " was used\n";
        } else {
            throw std::runtime_error("The primes file was too short");
        }
        primes.close();
    } else {
        try {
            primes.close();
        } catch (...) {}
        throw std::runtime_error("Could not open " + std::string(primes_source));
    }

    std::ifstream input((std::string(seeds_source)));
    std::string property;
    if (input.is_open()) {
        while (!input.eof()) {
            input >> property;
            if (property == "RANDOMSEED") { input >> m_l1 >> m_l2 >> m_l3 >> m_l4; }
        }
        input.close();
    } else {
        try {
            input.close();
        } catch (...) {}
        throw std::runtime_error("Could not open " + std::string(seeds_source));
    }
}

/**
 * Generates a pseudorandom unsigned integer in the interval [0, (2^48)-1]
 * @return
 */
ARandom::result_type ARandom ::operator()() {
    result_type i1 = m_l1 * m_m4 + m_l2 * m_m3 + m_l3 * m_m2 + m_l4 * m_m1 + m_p1;
    result_type i2 = m_l2 * m_m4 + m_l3 * m_m3 + m_l4 * m_m2 + m_p2;
    result_type i3 = m_l3 * m_m4 + m_l4 * m_m3 + m_p3;
    result_type i4 = m_l4 * m_m4 + m_p4;
    m_l4 = i4 % 4096UL;
    i3 = i3 + i4 / 4096UL;
    m_l3 = i3 % 4096UL;
    i2 = i2 + i3 / 4096UL;
    m_l2 = i2 % 4096UL;
    m_l1 = (i1 + i2 / 4096UL) % 4096UL;
    return b4096tob10<result_type>(m_l1, m_l2, m_l3, m_l4);
}


[[maybe_unused]] void ARandom ::SaveSeed() const {
    std::ofstream WriteSeed;
    WriteSeed.open("seed.out");
    if (WriteSeed.is_open()) {
        WriteSeed << m_l1 << " " << m_l2 << " " << m_l3 << " " << m_l4 << std::endl;
    } else {
        try {
            WriteSeed.close();
        } catch (...) {}
        throw std::runtime_error("PROBLEM: Unable to open seed.out");
    }
    WriteSeed.close();
}

#ifdef __apple_build_version__
[[maybe_unused]] double ARandom ::Gauss(double mean, double sigma) {
    double s = Rannyu();
    double t = Rannyu();
    double x = sqrt(-2 * log(1 - s)) * cos(2 * M_PI * t);
    return mean + x * sigma;
}

double ARandom ::Rannyu(double min, double max) { return min + (max - min) * Rannyu(); }

double ARandom ::Rannyu() {
    int i1 = m_l1 * m_m4 + m_l2 * m_m3 + m_l3 * m_m2 + m_l4 * m_m1 + m_p1;
    int i2 = m_l2 * m_m4 + m_l3 * m_m3 + m_l4 * m_m2 + m_p2;
    int i3 = m_l3 * m_m4 + m_l4 * m_m3 + m_p3;
    int i4 = m_l4 * m_m4 + m_p4;
    m_l4 = i4 % 4096;
    i3 = i3 + i4 / 4096;
    m_l3 = i3 % 4096;
    i2 = i2 + i3 / 4096;
    m_l2 = i2 % 4096;
    m_l1 = (i1 + i2 / 4096) % 4096;
    return twom12 * (double(m_l1) +
                     twom12 * (double(m_l2) + twom12 * (double(m_l3) + twom12 * double(m_l4))));
}
#endif

[[maybe_unused]] void ARandom ::SetRandom(result_type const *s, result_type p1, result_type p2) {
    m_l1 = s[0];
    m_l2 = s[1];
    m_l3 = s[2];
    m_l4 = s[3];
    m_p3 = p1;
    m_p4 = p2;
}


void ARandom::seed(ARandom::result_type s) { b10tob4096(s, &m_l1, &m_l2, &m_l3, &m_l4); }

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
