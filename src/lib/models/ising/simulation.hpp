//
// Created by Davide Nicoli on 04/05/22.
//

#ifndef ESERCIZI_LSN_LIBS_MODELS_ISING_SIMULATOR_HPP
#define ESERCIZI_LSN_LIBS_MODELS_ISING_SIMULATOR_HPP

#include <iostream>
#include <memory>
#include <string_view>

#include "ariel_random/ariel_random.hpp"
#include "ising.hpp"
#include "samplers/MCMC/gibbs.hpp"
#include "samplers/MCMC/metropolis.hpp"
#include "structs.hpp"

using samplers::mcmc::Gibbs;
using samplers::mcmc::Metropolis;

namespace ising {


    /*
template <typename var_space, typename p_space> class Simulator {
  using IsingModel = Ising<1, var_space>; // NOLINT(clion-misra-cpp2008-5-0-4)
  using IsingState = typename IsingModel::State;
  using Spin = typename IsingState::value_type;
  using MetropolisSampler = Metropolis<var_space, IsingState, UniformSpinFlip<p_space>>;
  using GibbsSampler = Gibbs<var_space, IsingState1D>;

private:
  StateVariables<var_space> m_thermo_state;
  const Parameters<var_space> m_parameters;
  const SimulationConfig m_simconf;
  IsingModel m_model;
  using SamplerPtr = std::unique_ptr<MCMCSampler<var_space, IsingState>>;
  SamplerPtr m_sampler;

public:
  Simulator(const std::string_view &config, ARandom &rng)
      : m_thermo_state(config), m_parameters(config), m_simconf(config),
        m_model(Tinf<1, var_space>(m_parameters.n_spins(), rng)) {
    if (m_simconf.is_metro()) {
      m_sampler = std::make_unique<MetropolisSampler>(
          m_model.m_state,
          [&](const IsingState &state) {
            return IsingModel::log_probability(state, m_parameters.J(),
                                               m_parameters.h(),
                                               1 / m_thermo_state.temperature());
          },
          m_parameters.n_spins(), true);
      std::cout << "Using metro" << '\n';
    } else {
      m_sampler = std::make_unique<GibbsSampler>(
          m_model.m_state,
          [&](const IsingState &state) {
            return IsingModel::log_probability(state, m_parameters.J(),
                                               m_parameters.h(),
                                               1 / m_thermo_state.temperature());
          });
    }
  }

  Simulator(Simulator &) = delete;
  Simulator(const Simulator &) = delete;

  auto &operator=(Simulator &) = delete;
  auto &operator=(const Simulator &) = delete;

  void step(ARandom &rng) {
    m_model.m_state = m_sampler->step(rng).second;
    m_thermo_state.internalEnergy = m_model.hamiltonian(m_parameters.J(), m_parameters.h());
    m_thermo_state.magnetization = std::reduce(m_model.state.cbegin(), m_model.state.cend(), spins::reduce, var_space{0});
    m_thermo_state.heatCapacity = 0;
  }
  const IsingState1D & state() const { return m_model.m_state; }
  const StateVariables<var_space>& thermo_state() const {return m_thermo_state;}
  const Parameters<var_space> & parameters() const { return m_parameters; }
};*/

}// namespace ising

#endif// ESERCIZI_LSN_LIBS_MODELS_ISING_SIMULATOR_HPP
