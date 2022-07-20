//
// Created by Davide Nicoli on 19/07/22.
//

#ifndef ESERCIZI_LSN_STRUCTURES_HPP
#define ESERCIZI_LSN_STRUCTURES_HPP

#include <array>
#include <cstddef>
#include <fstream>
#include <system_error>
#include <valarray>

using std::slice;

template<typename field>
struct Vectors {
    typedef std::valarray<field> v_type;

    explicit Vectors(size_t size) : e_i(size), e_j(size), e_k(size) {}
    // Vectors(size_t size) : m_size{size}, m_data(3 * size) {}

    // clang-format off
    Vectors(v_type &&v_i, v_type &&v_j, v_type &&v_k)
        : e_i{std::forward<v_type>(v_i)},
          e_j{std::forward<v_type>(v_j)},
          e_k{std::forward<v_type>(v_k)} {}
    // clang-format on

    /*Vectors(v_type &&v_i, v_type &&v_j, v_type &&v_k) : Vectors(v_i.size()) {
        m_data[slice(0, m_size, 1)] = std::forward<v_type>(v_i);
        m_data[slice(m_size, 2 * m_size, 1)] = std::forward<v_type>(v_j);
        m_data[slice(2 * m_size, 3 * m_size, 1)] = std::forward<v_type>(v_k);
    }*/

    Vectors(size_t size, std::string_view configuration_path) : Vectors(size) {
        std::ifstream configuration((std::string(configuration_path)));
        if (!configuration.is_open())
            throw std::runtime_error("Could not open " + std::string(configuration_path));

        for (size_t i = 0UL; i < size; i++) configuration >> e_i[i] >> e_j[i] >> e_k[i];
        configuration.close();
    }

    void save_configuration(std::string_view configuration_path) const {
        std::ofstream configuration((std::string(configuration_path)));
        if (!configuration.is_open())
            throw std::runtime_error("Could not open " + std::string(configuration_path));

        for (int i = 0; i < e_i.size(); ++i) {
            configuration << e_i[i] << "   " << e_j[i] << "   " << e_k[i] << std::endl;
        }
        configuration.close();
    }

    template<typename Transform>
    void save_configuration(std::string_view configuration_path, Transform f) const {
        std::ofstream configuration((std::string(configuration_path)));
        if (!configuration.is_open())
            throw std::runtime_error("Could not open " + std::string(configuration_path));

        for (int i = 0; i < e_i.size(); ++i) {
            configuration << f(e_i[i]) << "   " << f(e_j[i]) << "   " << f(e_k[i]) << std::endl;
        }
        configuration.close();
    }

    friend Vectors operator-(const Vectors &lhs, const Vectors &rhs) noexcept {
        return {lhs.e_i - rhs.e_i, lhs.e_j - rhs.e_j, lhs.e_k - rhs.e_k};
    }

    Vectors &operator-=(const std::array<field, 3> &oth) noexcept {
        e_i -= oth[0];
        e_j -= oth[1];
        e_k -= oth[2];
        return *this;
    }
    Vectors<field> operator-(const std::array<field, 3> &oth) const noexcept {
        return {e_i - oth[0], e_j - oth[1], e_k - oth[2]};
    }
    Vectors &operator-=(const Vectors<field> &oth) noexcept {
        e_i -= oth.e_i;
        e_j -= oth.e_j;
        e_k -= oth.e_k;
        return *this;
    }


    Vectors &operator*=(field x) noexcept {
        e_i *= x;
        e_j *= x;
        e_k *= x;
        return *this;
    }

    std::array<field, 3> operator[](size_t i) const { return {e_i[i], e_j[i], e_k[i]}; }
    template<size_t I>
    field at(size_t i) {
        static_assert(I < 3);
        if constexpr (I == 0) {
            return e_i[i];
        } else if constexpr (I == 1) {
            return e_j[i];
        } else if constexpr (I == 2) {
            return e_k[i];
        }
    }

    std::array<field, 3> mean() const noexcept {
        return {e_i.sum() / field(e_i.size()), e_j.sum() / field(e_j.size()),
                e_k.sum() / field(e_k.size())};
    }

    field full_norm2() const noexcept { return (e_i * e_i + e_j * e_j + e_k * e_k).sum(); }
    v_type axial_norm2() const noexcept { return e_i * e_i + e_j * e_j + e_k * e_k; };

    template<typename Transform>
    Vectors &transform(Transform f) {
        std::for_each(std::begin(e_i), std::end(e_i), f);
        std::for_each(std::begin(e_j), std::end(e_j), f);
        std::for_each(std::begin(e_k), std::end(e_k), f);
        return *this;
    }

    template<typename Transform>
    Vectors<field> apply(Transform f) {
        return Vectors<field>(e_i.apply(f), e_j.apply(f), e_k.apply(f));
    }

    v_type e_i, e_j, e_k;
    // size_t m_size;
    // v_type m_data;
};

template<typename field>
struct ThermoSettings {
    ThermoSettings(field T, field rho)
        : temperature{T}, density{rho}, beta{field(1) / temperature} {}

    explicit ThermoSettings(std::string_view settings_path) {
        std::ifstream settings((std::string(settings_path)));
        if (!settings.is_open())
            throw std::runtime_error("Could not open " + std::string(settings_path));
        bool _bucket, _resume;
        field T, rho, r, dt;
        size_t _n_part;
        settings >> _bucket >> _resume >> T >> _n_part >> rho;
        settings.close();
        ThermoSettings<field>(T, rho);
    }

    field temperature, density, beta;
};

template<typename field>
struct SimulationSettings {
    SimulationSettings(size_t n_part, size_t n_b, size_t n_s, field r, field deltat, field density)
        : n_particles{n_part}, n_blocks{n_b}, n_steps{n_s}, cutoff{r}, dt{deltat},
          volume{field(n_particles) / density}, box_edge{std::cbrt(volume)} {}

    explicit SimulationSettings(std::string_view settings_path) {
        std::ifstream settings((std::string(settings_path)));
        if (!settings.is_open())
            throw std::runtime_error("Could not open " + std::string(settings_path));
        bool _bucket, _resume;
        field T, rho, r, deltat;
        size_t n_part, n_b, n_s;
        settings >> _bucket >> _resume >> T >> n_part >> rho >> r >> deltat >> n_b >> n_s;
        settings.close();
        SimulationSettings<field>(n_part, n_b, n_s, r, deltat, rho);
    }

    size_t n_particles, n_blocks, n_steps;
    field cutoff, dt;
    field volume, box_edge;
};

struct ThermoVariables {};

#endif//ESERCIZI_LSN_STRUCTURES_HPP
