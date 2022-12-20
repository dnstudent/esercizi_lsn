//
// Created by Davide Nicoli on 03/05/22.
//

#ifndef ESERCIZI_LSN_ISING_STRUCTS_HPP
#define ESERCIZI_LSN_ISING_STRUCTS_HPP

#include <exception>
#include <fstream>
#include <string>
#include <string_view>
#include <valarray>
#include <vector>

/*template<typename real>
using Vector5 = std::valarray<real>;*/

namespace ising {
    namespace spins {
        typedef bool BinarySpin;

        constexpr bool up = true;
        constexpr int up_int = 1;
        constexpr bool down = false;
        constexpr int down_int = -1;

        /**
         * Spin multiplication. (u * d) = (d * u) = -1; (u * u) = (d * d) = +1;
         * @tparam integer Return type.
         * @param a First spin
         * @param b Second spin
         */
        template<typename integer>
        integer multiplies(const BinarySpin a, const BinarySpin b) {
            return (a ^ b) ? -1 : 1;
        }

        /**
         * Spin addition. (u + d) = (d + u) = 0; (u + u) = +2; (d + d) = -2;
         * @tparam integer Return type.
         * @param a First spin.
         * @param b Second spin.
         */
        template<typename integer>
        integer plus(const bool a, const bool b) {
            integer result;
            if (a ^ b) {
                result = 0;
            } else if (a) {
                result = 2;
            } else {
                result = -2;
            }
            return result;
        }

        /**
         * The sum reductor: it can be applied with an accumulating function (e.g. std::accumulate) to sum over a vector of spins.
         * WATCH OUT: it is not commutative nor associative! It cannot be used with std::reduce.
         * @tparam integer Return type.
         * @param acc Accumulator.
         * @param a Spin.
         * @return if the spin is up, the accumulator incremented, otherwise reduced.
         */
        template<typename integer>
        integer reductor(integer acc, const BinarySpin a) {
            if (a) { return acc + 1; }
            return acc - 1;
        }
    };// namespace spins
}// namespace ising

#endif// ESERCIZI_LSN_ISING_STRUCTS_HPP
