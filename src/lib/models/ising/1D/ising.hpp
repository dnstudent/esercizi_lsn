//
// Created by Davide Nicoli on 29/08/22.
//

#ifndef ESERCIZI_LSN_ISING_1D_HPP
#define ESERCIZI_LSN_ISING_1D_HPP

#include <cassert>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <random>


#include "../ising.hpp"
#include "distributions/uniform_int.hpp"
#include "structs.hpp"
#include "variables.hpp"

namespace fs = std::filesystem;

namespace ising {
    template<typename var_space>
    class Ising<var_space, 1> {

    public:
        typedef var_space VarSpace;
        typedef IsingState1D StateSpace;

        template<class URBG>
        Ising(size_t n_spins, URBG &rng, VarSpace J, VarSpace h, VarSpace T)
            : m_state(n_spins), m_n_spins(n_spins), m_choosespin(0, n_spins - 1), m_J(J), m_h(h),
              m_beta(1 / T) {
            assert(n_spins > 2);
            std::uniform_real_distribution<double> unif(0, 1);
            std::generate(m_state.begin(), m_state.end(),
                          [&]() -> spins::BinarySpin { return unif(rng) < 0.5; });
        }

        Ising(size_t n_spins, spins::BinarySpin spin_value, VarSpace J, VarSpace h, VarSpace T)
            : m_state(n_spins, spin_value), m_n_spins(n_spins), m_choosespin(0, n_spins - 1),
              m_J(J), m_h(h), m_beta(1 / T) {
            assert(n_spins > 2);
        }

        Ising(std::initializer_list<spins::BinarySpin> state, VarSpace J, VarSpace h, VarSpace T)
            : m_state(state), m_n_spins(m_state.size()), m_choosespin(0, m_n_spins - 1), m_J(J),
              m_h(h), m_beta(1 / T) {
            assert(m_n_spins > 2);
        }

        Ising(const fs::path &state_path, VarSpace J, VarSpace h, VarSpace T)
            : m_state(read_state(state_path)), m_n_spins(m_state.size()),
              m_choosespin(0, m_n_spins - 1), m_J(J), m_h(h), m_beta(1 / T) {
            assert(m_n_spins > 2);
        }

        spins::BinarySpin operator[](int64_t k) const {
            while (k < 0) k = m_n_spins_signed + k;
            return m_state[unsigned(k % m_n_spins_signed)];
        }

        void set(size_t k, spins::BinarySpin value) { m_state[k] = value; }

        inline void flip(size_t k) noexcept { m_state[k] = !m_state[k]; }

        template<class URBG>
        inline size_t sample(URBG &rng) noexcept {
            return m_choosespin(rng);
        }

        inline void evolve(size_t spin) noexcept { flip(spin); }

        [[nodiscard]] StateSpace const &state() const noexcept { return m_state; }


        // ∆E when flipping spin k
        inline var_space flip_dE(size_t k) const noexcept {
            const auto s_k = m_state[k] ? spins::up_int : spins::down_int;
            return 2 * s_k *
                   (m_J * spins::plus<var_space>((*this)[int64_t(k) - 1], (*this)[int64_t(k) + 1]) +
                    m_h);
        }


        var_space energy() const noexcept {
            const auto s0 = m_state.cbegin();
            const auto s1 = std::next(s0);
            const auto sN = std::prev(m_state.cend());

            const auto _h1 =
                    var_space(std::transform_reduce(s0, sN, s1, spins::multiplies<int>(*sN, *s0),
                                                    std::plus<>(), spins::multiplies<int>));
            if (m_h == var_space(0)) return -m_J * _h1;
            return -m_J * _h1 - m_h * D1::spin_sum(*this);
        }


        inline var_space logp() const noexcept { return -m_beta * energy(); }

        template<typename prob_space>
        inline prob_space flip_logp(size_t candidate) const noexcept {
            return -m_beta * flip_dE(candidate);
        }


        void save_state(const fs::path &output_path) const {
            if (!fs::exists(output_path.parent_path()))
                fs::create_directories(output_path.parent_path());
            std::ofstream file(output_path);
            if (!file.is_open()) throw std::runtime_error("Could not open " + output_path.string());
            file << "State\n";
            for (size_t i = 0; i + 1 < m_n_spins; i++) { file << m_state[i] << '\n'; }
            file << m_state.back();
            file.close();
        }

        StateSpace read_state(const fs::path &state_path) {
            std::ifstream file(state_path);
            if (!file.is_open())
                throw std::runtime_error(state_path.string() + " could not be opened.");
            std::string _bucket;
            file >> _bucket;
            StateSpace state{};
            bool s;
            while (!file.eof()) {
                file >> s;
                state.push_back(s);
            }
            file.close();
            return state;
        }

        VarSpace h() const { return m_h; }

        VarSpace beta() const noexcept { return m_beta; }

        [[nodiscard]] size_t n_spins() const noexcept { return m_n_spins; }

    private:
        StateSpace m_state;
        size_t m_n_spins;
        int64_t m_n_spins_signed{signed(m_n_spins)};
        distributions::uniform_int<size_t> m_choosespin;
        const VarSpace m_J, m_h, m_beta;
    };

    namespace D1 {
        template<typename var_space, class URBG>
        Ising<var_space, 1> T0(size_t n_spins, var_space J, var_space h, var_space T, URBG &rng) {
            std::bernoulli_distribution coin(0.5);
            return {n_spins, coin(rng), J, h, T};
        }

        template<typename var_space, class URBG>
        Ising<var_space, 1> Tinf(size_t n_spins, var_space J, var_space h, var_space T, URBG &rng) {
            return {n_spins, rng, J, h, T};
        }
    }// namespace D1
}// namespace ising

#endif//ESERCIZI_LSN_ISING_1D_HPP
