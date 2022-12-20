//
// Created by Davide Nicoli on 20/10/22.
//

#ifndef ESERCIZI_LSN_MC_INTEGRATOR_HPP
#define ESERCIZI_LSN_MC_INTEGRATOR_HPP

#include <type_traits>
#include <utility>
#include <vector>

#include "estimators/estimators.hpp"
#include "estimators/mean.hpp"

using namespace estimators;

namespace mc_integrator {
    /**
     * Monte Carlo integrator
     * @tparam MCSampler Monte Carlo sampler (e.g. Metropolis)
     */
    template<typename field, class MCSampler>
    class Integrator {
        using X = typename MCSampler::StateSpace;

    public:
        /**
         * Constructor. The sampler must sample by the distributional part of the integral.
         */
        explicit Integrator(MCSampler &&sampler) : m_sampler(std::forward<MCSampler>(sampler)) {}

        /**
         * Integrate f using the provided sampler.
         * @param f Function to integrate. Must be written so that the probability distribution used by the sampler is normalized.
         * @param n_samples Number of sampled points.
         * @param rng Random number generator.
         * @return Estimate.
         */
        template<class Integrand, class URBG>
        auto operator()(Integrand f, size_t n_samples, URBG &rng) {
            using f_out = std::invoke_result_t<Integrand, const X &>;
            std::vector<X> xs(n_samples);
            std::vector<f_out> ys(n_samples);

            m_sampler.sample(xs.begin(), xs.end(), rng);
            std::transform(xs.cbegin(), xs.cend(), ys.begin(), f);
            return utils::average<f_out>(ys.cbegin(), ys.cend());
        }

        /**
         * Integrate f using the provided sampler, returning both the estimate and the uncertainty.
         * @param f Function to integrate. Must be written so that the probability distribution used by the sampler is normalized.
         * @param n_blocks Number of data blocks.
         * @param block_size Number of MC steps in each block.
         * @param rng Random number generator.
         * @return Pair estimate, statistical uncertainty.
         */
        template<class Integrand, class URBG>
        auto operator()(Integrand f, size_t n_blocks, size_t block_size, URBG &rng) {
            using f_out = std::invoke_result_t<Integrand, const X &>;
            std::vector<X> xs(block_size);
            std::vector<f_out> ys(block_size);
            // Performing the necessary steps, throwing away the results
            for (size_t block = 0; block + 1 < n_blocks; block++) {
                m_sampler.sample(xs.begin(), xs.end(), rng);
                std::transform(xs.cbegin(), xs.cend(), ys.begin(), f);
                m_estimator(ys.cbegin(), ys.cend());
            }
            m_sampler.sample(xs.begin(), xs.end(), rng);
            std::transform(xs.cbegin(), xs.cend(), ys.begin(), f);
            return m_estimator(ys.cbegin(), ys.cend());
        }

        template<class Integrand, typename EstimateOut, typename UncertOut, typename SpaceOut,
                 class URBG>
        void operator()(Integrand f, EstimateOut first_estimate, EstimateOut last_estimate,
                        UncertOut first_uncert, SpaceOut first_x, SpaceOut last_x, URBG &rng) {
            const auto n_blocks = static_cast<size_t>(std::distance(first_estimate, last_estimate));
            if (static_cast<size_t>(std::distance(first_x, last_x)) % n_blocks != 0)
                throw std::runtime_error("Not divisible");
            const auto block_size = static_cast<size_t>(std::distance(first_x, last_x)) / n_blocks;
            using f_out = std::invoke_result_t<Integrand, const X &>;
            std::vector<f_out> ys(block_size);
            SpaceOut next_x;
            while (first_x != last_x) {
                next_x = utils::snext(first_x, block_size);
                m_sampler.sample(first_x, next_x, rng);
                std::transform(first_x, next_x, ys.begin(), f);
                std::tie(*first_estimate++, *first_uncert++) = m_estimator(ys.cbegin(), ys.cend());
                first_x = next_x;
            }
        }

        /**
         * Integrate f to a provided statistical uncertainty.
         * @param f Function to integrate. Must be written so that the probability distribution used by the sampler is normalized.
         * @param statistical_uncertainty Statistical uncertainty to integrate to.
         * @param block_size Number of MC steps in each block.
         * @param rng Random number generator.
         * @return Pair (Pair (estimate, statistical uncertainty), number of steps).
         */
        template<class Integrand, class URBG>
        std::pair<std::tuple<field, field>, size_t>
        integrate_to(Integrand f, field statistical_uncertainty, size_t block_size, URBG &rng) {
            using f_out = std::invoke_result_t<Integrand, const X &>;
            std::vector<X> xs(block_size);
            std::vector<f_out> ys(block_size);
            m_sampler.sample(xs.begin(), xs.end(), rng);
            std::transform(xs.cbegin(), xs.cend(), ys.begin(), f);
            auto result = m_estimator(ys.cbegin(), ys.cend());
            size_t n_blocks = 0;
            do {
                m_sampler.sample(xs.begin(), xs.end(), rng);
                std::transform(xs.cbegin(), xs.cend(), ys.begin(), f);
                result = m_estimator(ys.cbegin(), ys.cend());
                n_blocks++;
            } while (std::get<1>(result) > statistical_uncertainty);
            return std::make_pair(result, n_blocks);
        }

    private:
        MCSampler m_sampler;
        ProgAvg<field> m_estimator{};
    };
}// namespace mc_integrator

#endif//ESERCIZI_LSN_MC_INTEGRATOR_HPP
