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
 * @param seeds_source Path to transform file containing the seed. Must be on a line starting with "RANDOMSEED" followed by four base-2^12 digits separated by spaces.
 * @param primes_source Path to a file containing primes. Each line of the file must list the two base-2^12 digits representing a prime < 2^24, separated by transform single space.
 * @param primes_line The line of the primes file to use, starting from 1. If 0 use the default prime (11848219)
 */
ARandom::ARandom(const std::string_view &seeds_source, const std::string_view &primes_source,
                 size_t primes_line) {
    read_primes(primes_source, primes_line, m_p3, m_p4);
    read_seeds(seeds_source, m_l1, m_l2, m_l3, m_l4);
}

/**
 * Generates transform pseudorandom unsigned integer in the interval [0, (2^48)-1]
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


void ARandom ::SaveSeed(std::string_view path) const {
    std::ofstream WriteSeed((std::string(path)));
    if (WriteSeed.is_open()) {
        WriteSeed << "RANDOMSEED"
                  << " " << m_l1 << " " << m_l2 << " " << m_l3 << " " << m_l4 << std::endl;
        WriteSeed.close();
    } else {
        try {
            WriteSeed.close();
        } catch (...) {}
        throw std::runtime_error("PROBLEM: Unable to open seed.out");
    }
}

[[maybe_unused]] double ARandom ::Gauss(double mean, double sigma) {
    double s = Rannyu();
    double t = Rannyu();
    double x = sqrt(-2 * log(1 - s)) * cos(2 * M_PI * t);
    return mean + x * sigma;
}

double ARandom ::Rannyu(double min, double max) { return min + (max - min) * Rannyu(); }

double ARandom ::Rannyu() {
    result_type i1 = m_l1 * m_m4 + m_l2 * m_m3 + m_l3 * m_m2 + m_l4 * m_m1 + m_p1;
    result_type i2 = m_l2 * m_m4 + m_l3 * m_m3 + m_l4 * m_m2 + m_p2;
    result_type i3 = m_l3 * m_m4 + m_l4 * m_m3 + m_p3;
    result_type i4 = m_l4 * m_m4 + m_p4;
    m_l4 = i4 % 4096;
    i3 = i3 + i4 / 4096;
    m_l3 = i3 % 4096;
    i2 = i2 + i3 / 4096;
    m_l2 = i2 % 4096;
    m_l1 = (i1 + i2 / 4096) % 4096;
    return twom12 * (double(m_l1) +
                     twom12 * (double(m_l2) + twom12 * (double(m_l3) + twom12 * double(m_l4))));
}

[[maybe_unused]] void ARandom ::SetRandom(result_type const *s, result_type p1, result_type p2) {
    m_l1 = s[0];
    m_l2 = s[1];
    m_l3 = s[2];
    m_l4 = s[3];
    m_p3 = p1;
    m_p4 = p2;
}


void ARandom::seed(ARandom::result_type s) { b10tob4096(s, &m_l1, &m_l2, &m_l3, &m_l4); }
void ARandom::save_seed(std::string_view path,
                        std::ios_base::openmode mode = std::ios_base::out) const {
    size_t s = b4096tob10(m_l1, m_l2, m_l3, m_l4);
    std::ofstream seed_file((std::string(path)), mode);
    if (!seed_file.is_open()) { throw std::runtime_error("Could not open " + std::string(path)); }
    seed_file << s << std::endl;
    seed_file.close();
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
