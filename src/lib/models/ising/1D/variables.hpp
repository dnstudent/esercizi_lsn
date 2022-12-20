//
// Created by Davide Nicoli on 30/08/22.
//

#ifndef ESERCIZI_LSN_ISING_1D_VARIABLES_HPP
#define ESERCIZI_LSN_ISING_1D_VARIABLES_HPP

#include <algorithm>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <vector>

#include <rapidcsv.h>

#include "../structs.hpp"
#include "estimators/mean.hpp"
#include "estimators/variance.hpp"
#include "ising.hpp"

namespace csv = rapidcsv;

namespace ising::D1 {

    /**
     * Struct for the proxy variables: H (the energy), sum_s (the sum of spins) and sum_s2 (the sum of the product of spins)
     * @tparam var_space The numeric field to use for H.
     */
    template<typename var_space>
    struct block_proxy_vars {
        explicit block_proxy_vars(size_t block_size)
            : h(block_size), sum_s(block_size), sum_s2(block_size) {}
        explicit block_proxy_vars(const fs::path &data) {
            if (!fs::exists(data))
                throw std::runtime_error(data.string() +
                                         " does not exist. Cwd: " + fs::current_path().string());
            csv::Document table(data);
            h = table.GetColumn<var_space>("H");
            sum_s = table.GetColumn<int64_t>("Sum_s");
            sum_s2 = table.GetColumn<int64_t>("Sum_s2");
        }
        std::vector<var_space> h;
        std::vector<int64_t> sum_s, sum_s2;

        void save_data(const fs::path &output_path) const {

            csv::Document table;
            table.InsertColumn(0, h, "H");
            table.InsertColumn(1, sum_s, "Sum_s");
            table.InsertColumn(2, sum_s2, "Sum_s2");
            table.RemoveColumn(3);

            if (!fs::exists(output_path.parent_path()))
                fs::create_directories(output_path.parent_path());
            table.Save(output_path);
        }
    };

    /**
     * Computes the sum of the spins of a 1D Ising model
     * @param model Ising model.
     * @return Spins sum: sum_s.
     */
    template<typename var_space>
    inline int spin_sum(const Ising<var_space, 1> &model) noexcept {
        return std::accumulate(model.state().cbegin(), model.state().cend(), int(0),
                               spins::reductor<int>);
    }

    /**
     * Computes the sum of the product of each spin pair in a 1D Ising model
     * @param model 1D Ising model.
     * @return The aforementioned sum: sum_s2.
     */
    template<typename var_space>
    inline int64_t spin_sum2(const Ising<var_space, 1> &model) noexcept {
        const auto ss = spin_sum(model);
        return ss * ss;
    }

    /**
     * Estimator of the internal energy per spin of a 1D Ising model
     * @tparam mean_estimator Which object should be used as an estimator of the mean.
     */
    template<typename var_space, class mean_estimator>
    class InternalEnergy {
    public:
        typedef std::pair<var_space, var_space> return_type;
        explicit InternalEnergy(Ising<var_space, 1> &ising) : n_spins(var_space(ising.n_spins())) {
            assert(ising.h() == 0);
        }

        /**
         * Computes an estimate of <H> / N using the mean_estimator on proxy variables.
         * @param block_data Proxy variables sample.
         * @return Estimate of the internal energy per spin and its uncertainty
         */
        inline constexpr return_type operator()(const block_proxy_vars<var_space> &block_data) {
            const auto [estimate, error] =
                    m_mean_estimator(block_data.h.cbegin(), block_data.h.cend());
            return {estimate / n_spins, error / n_spins};
        }

        [[nodiscard]] std::string name() const { return "u"; }

    private:
        var_space n_spins;
        mean_estimator m_mean_estimator{};
    };

    /**
     * Estimator of the heat capacity per spin of a 1D Ising model.
     * @tparam variance_estimator Estimator used as an estimator of the variance.
     */
    template<typename var_space, class variance_estimator>
    class HeatCapacity {
    public:
        typedef std::pair<var_space, var_space> return_type;
        explicit HeatCapacity(Ising<var_space, 1> &ising)
            : m_coeff(std::pow(ising.beta(), 2) / var_space(ising.n_spins())),
              m_coeff2(std::pow(ising.beta(), 4) / var_space(ising.n_spins())) {
            assert(ising.h() == 0);
        }

        /**
         * Computes an estimate of b^2*<<H>> / N using the variance_estimator on proxy variables.
         * @param block_data Proxy variables sample.
         * @return Estimate of the heat capacity per spin and its uncertainty
         */
        inline constexpr return_type operator()(const block_proxy_vars<var_space> &block_data) {
            const auto [mean, error] =
                    m_variance_estimator(block_data.h.cbegin(), block_data.h.cend());
            return {m_coeff * mean, m_coeff2 * error};
        }

        [[nodiscard]] std::string name() const { return "c"; }

    private:
        var_space m_coeff, m_coeff2;
        variance_estimator m_variance_estimator{};
    };

    /**
     * Estimator of the magnetic susceptivity of a 1D Ising model.
     * @tparam mean_estimator Estimator used as an estimator of the mean.
     */
    template<typename var_space, class mean_estimator>
    class MagneticSusceptivity {
    public:
        typedef std::pair<var_space, var_space> return_type;
        explicit MagneticSusceptivity(Ising<var_space, 1> &ising)
            : m_coeff(ising.beta() / var_space(ising.n_spins())), m_coeff2(m_coeff * ising.beta()) {
            assert(ising.h() == 0);
        }

        /**
         * Computes an estimate of b*<sum_s2> / N using the mean_estimator on proxy variables.
         * @param block_data Proxy variables sample.
         * @return Estimate of the magnetic susceptivity and its uncertainty
         */
        inline constexpr return_type operator()(const block_proxy_vars<var_space> &block_data) {
            const auto [mean, error] =
                    m_mean_estimator(block_data.sum_s2.cbegin(), block_data.sum_s2.cend());
            return {m_coeff * var_space(mean), m_coeff2 * var_space(error)};
        }
        [[nodiscard]] std::string name() const { return "X"; }

    private:
        var_space m_coeff, m_coeff2;
        mean_estimator m_mean_estimator{};
    };

    /**
     * Estimator of the magnetization per spin of a 1D Ising model.
     * @tparam mean_estimator Estimator used as an estimator of the mean.
     */
    template<typename var_space, class mean_estimator>
    class Magnetization {
    public:
        typedef std::pair<var_space, var_space> return_type;
        explicit Magnetization(Ising<var_space, 1> &ising) : n_spins(var_space(ising.n_spins())) {}

        /**
         * Computes an estimate of <sum_s> / N using the mean_estimator on proxy variables.
         * @param block_data Proxy variables sample.
         * @return Estimate of the magnetization per spin and its uncertainty
         */
        inline constexpr return_type operator()(const block_proxy_vars<var_space> &block_data) {
            const auto [estimate, error] =
                    m_mean_estimator(block_data.sum_s.cbegin(), block_data.sum_s.cend());
            return {var_space(estimate) / n_spins, var_space(error) / n_spins};
        }
        [[nodiscard]] std::string name() const { return "m"; }

    private:
        var_space n_spins;
        mean_estimator m_mean_estimator{};
    };

    /**
     * Computes and stores estimates of the provided variables using proxy samples
     * @tparam I Internal use.
     * @tparam Variables Tuple-like object which members are variables of the 1D Ising model.
     * @tparam Outs Tuple-like object which members are pair of vectors that will store respectively estimates and uncertainties. A proper object, for example, is array<pair<vector<var_space>, vector<var_space>>, tuple_size_v<Variables>>
     */
    template<size_t I = 0, typename var_space, class Variables, class Outs>
    void compute_store_estimates(block_proxy_vars<var_space> &block_data, Variables &variables,
                                 Outs &outs) {
        // This is a compile time recursive function: it is necessary to iterate over a generic (though defined at compile time)
        // group of variables, which must be stored in a heterogeneus container (e.g. tuple).
        static_assert(std::tuple_size_v<Outs> == std::tuple_size_v<Variables>);
        if constexpr (I == std::tuple_size_v<Outs>) {
            return;
        } else {
            // Computation of the Ith variable estimate
            const auto [estimate, error] = std::get<I>(variables)(block_data);
            // Estimate and error are stored
            std::get<I>(outs).first.push_back(estimate);
            std::get<I>(outs).second.push_back(error);
            // recursion: it is not possible to iterate over a tuple differently (e.g. with a for loop)
            compute_store_estimates<I + 1>(block_data, variables, outs);
        }
    }
}// namespace ising::D1

#endif//ESERCIZI_LSN_ISING_1D_VARIABLES_HPP
