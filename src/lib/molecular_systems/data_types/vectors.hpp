//
// Created by Davide Nicoli on 28/09/22.
//

#ifndef ESERCIZI_LSN_MS_VECTORS_HPP
#define ESERCIZI_LSN_MS_VECTORS_HPP

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <utility>
#include <valarray>

namespace fs = std::filesystem;

/**
 * Object representing 3D coordinates as three lists
 * @tparam field The numeric field of the coordinates
 */
template<typename field>
struct Vectors {
    typedef std::valarray<field> v_type;

    explicit Vectors(size_t size) : e_i(size), e_j(size), e_k(size) {}

    // clang-format off
    Vectors(v_type &&v_i, v_type &&v_j, v_type &&v_k)
        : e_i{std::forward<v_type>(v_i)},
          e_j{std::forward<v_type>(v_j)},
          e_k{std::forward<v_type>(v_k)} {}
    // clang-format on

    /**
     * Initialize from transform configuration file consisting of three columns of numbers separated by spaces
     * @param size
     * @param configuration_path
     */
    Vectors(size_t size, const fs::path &configuration_path) : Vectors(size) {
        std::ifstream configuration(configuration_path);
        if (!configuration.is_open())
            throw std::runtime_error("Could not open " + configuration_path.string());

        for (size_t i = 0UL; i < size; i++) configuration >> e_i[i] >> e_j[i] >> e_k[i];
        configuration.close();
    }


    /**
     * Store the current state (transformed) as three columns of numbers separated by spaces
     * @param configuration_path The path where the configuration will be stored
     * @param f A transformation to be applied to the coordinates. Takes a scalar value and returns transform scalar value
     */
    template<typename Transform>
    void save_configuration(const fs::path &configuration_path, Transform f) const {
        std::ofstream configuration(configuration_path);
        if (!configuration.is_open())
            throw std::runtime_error("Could not open " + configuration_path.string());

        for (size_t i = 0; i < e_i.size(); ++i) {
            configuration << f(e_i[i]) << "   " << f(e_j[i]) << "   " << f(e_k[i]) << std::endl;
        }
        configuration.close();
    }

    /**
         * Store the current state (raw) as three columns of numbers separated by spaces
         * @param configuration_path The path where the configuration will be stored
         */
    inline void save_configuration(const fs::path &configuration_path) const {
        save_configuration(configuration_path, [](const auto x) { return x; });
    }
    /**
     * Store the current state (transformed) in the xyz format
     * @param configuration_path The path where the configuration will be stored
     * @param f A transformation to be applied to the coordinates. Takes a scalar value and returns transform scalar value
     */
    template<typename Transform>
    void save_xyz_configuration(const fs::path &configuration_path, Transform f) const {
        std::ofstream configuration(configuration_path);
        if (!configuration.is_open())
            throw std::runtime_error("Could not open " + configuration_path.string());
        configuration << e_i.size() << std::endl;
        configuration << "Comment" << std::endl;
        for (size_t i = 0; i < e_i.size(); ++i) {
            configuration << "LJ  " << f(e_i[i]) << "   " << f(e_j[i]) << "   " << f(e_k[i])
                          << std::endl;
        }
        configuration.close();
    }
    

    friend Vectors operator-(const Vectors &lhs, const Vectors &rhs) noexcept {
        return {lhs.e_i - rhs.e_i, lhs.e_j - rhs.e_j, lhs.e_k - rhs.e_k};
    }

    Vectors &operator=(const Vectors<field> &oth) {
        e_i.resize(oth.e_i.size());
        std::copy(std::begin(oth.e_i), std::end(oth.e_i), std::begin(e_i));
        e_j.resize(oth.e_j.size());
        std::copy(std::begin(oth.e_j), std::end(oth.e_j), std::begin(e_j));
        e_k.resize(oth.e_k.size());
        std::copy(std::begin(oth.e_k), std::end(oth.e_k), std::begin(e_k));
        return *this;
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

    friend bool operator==(const Vectors<field> &lhs, const Vectors<field> &rhs) {
        return std::equal(std::begin(lhs.e_i), std::end(lhs.e_i), std::begin(rhs.e_i)) &&
               std::equal(std::begin(lhs.e_j), std::end(lhs.e_j), std::begin(rhs.e_j)) &&
               std::equal(std::begin(lhs.e_k), std::end(lhs.e_k), std::begin(rhs.e_k));
    }


    Vectors &operator*=(field x) noexcept {
        e_i *= x;
        e_j *= x;
        e_k *= x;
        return *this;
    }

    /**
     * Retrieve a particle's coordinates
     * @param i
     * @return
     */
    std::array<field, 3> operator[](size_t i) const { return {e_i[i], e_j[i], e_k[i]}; }

    std::array<field, 3> mean() const noexcept {
        return {e_i.sum() / field(e_i.size()), e_j.sum() / field(e_j.size()),
                e_k.sum() / field(e_k.size())};
    }

    /**
     * Sum of all coordinates squared
     * @return
     */
    field full_norm2() const noexcept { return (e_i * e_i + e_j * e_j + e_k * e_k).sum(); }

    /**
     * Applies f to each coordinate in-place
     * @param f Transformation. Takes a scalar and returns transformed scalar
     */
    template<typename Transform>
    Vectors &apply(Transform f) {
        std::for_each(std::begin(e_i), std::end(e_i), f);
        std::for_each(std::begin(e_j), std::end(e_j), f);
        std::for_each(std::begin(e_k), std::end(e_k), f);
        return *this;
    }

    /**
     * Returns a vector of transformed coordinates
     * @param f Transformation. Takes a scalar and returns transformed scalar
     */
    template<typename Transform>
    Vectors<field> transform(Transform f) {
        return Vectors<field>(e_i.apply(f), e_j.apply(f), e_k.apply(f));
    }

    v_type e_i, e_j, e_k;
};


#endif//ESERCIZI_LSN_MS_VECTORS_HPP
