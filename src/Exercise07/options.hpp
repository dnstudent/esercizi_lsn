//
// Created by Davide Nicoli on 19/09/22.
//

#ifndef ESERCIZI_LSN_07_OPTS_HPP
#define ESERCIZI_LSN_07_OPTS_HPP

#include <filesystem>
#include <iostream>
#include <optional>

#include <cxxopts.hpp>

#include "utils.hpp"

namespace fs = std::filesystem;

namespace ex07 {
    /**
     * Struct storing program options for exercise 07.2
     * @tparam var_space The numeric field to use for the variables
     */
    struct Ex2Options {
        explicit Ex2Options(cxxopts::ParseResult &pr)
            : output_dir(pr["out"].as<fs::path>()), input_dir(pr["in"].as<fs::path>()),
              positions_path(input_dir / "positions"), velocities_path(input_dir / "velocities"),
              settings_path(pr["settings"].as<std::string>().empty()
                                    ? input_dir / "input"
                                    : fs::path(pr["settings"].as<std::string>())),
              output_positions(output_dir / "positions"),
              output_velocities(output_dir / "velocities"), output_settings(output_dir / "input"),
              save_every(pr["N"].as<size_t>()), resume(fs::exists(velocities_path)) {
            utils::require_existence(input_dir);
            utils::require_existence(positions_path);
            utils::require_existence(settings_path);
            if (!fs::exists(output_dir)) { fs::create_directories(output_dir); }
            if (!fs::exists(output_settings)) fs::copy_file(settings_path, output_settings);
            rng_seed_path = fs::exists(input_dir / "rng.seed")
                                    ? input_dir / "rng.seed"
                                    : fs::path(pr["s"].as<std::string>());
        }

        fs::path output_dir, input_dir, positions_path, velocities_path, settings_path,
                output_positions, output_velocities, output_settings, rng_seed_path;
        size_t save_every;
        bool resume;
    };

    enum Method { MC, MD };
    std::string tag(Method m) {
        if (m == MC) return "mc";
        else if (m == MD)
            return "md";
        else
            throw std::runtime_error("Not an option");
    }

    /**
     * Struct storing program options for exercise 07.4
     * @tparam var_space The numeric field to use for the variables
     */
    struct Ex4Options {
        explicit Ex4Options(cxxopts::ParseResult &pr)
            : input_dir({pr["in_mc"].as<fs::path>(), pr["in_md"].as<fs::path>()}),
              output_dir({pr["out"].as<fs::path>() / tag(MC), pr["out"].as<fs::path>() / tag(MD)}),
              n_bins(pr["n"].as<size_t>()), sample({pr["mc"].as<bool>(), pr["md"].as<bool>()}),
              warmup(pr["warmup"].as<bool>()) {
            for (auto m: {MC, MD}) {
                const auto input = pr[tag(m) + "_settings"].as<std::string>();
                input_settings[m] = input.empty() ? input_dir[m] / "input" : fs::path(input);
                if (!sample[m]) continue;
                utils::require_existence(input_settings[m]);
                utils::require_existence(input_positions[m]);
                if (!fs::exists(output_dir[m])) { fs::create_directories(output_dir[m]); }
                if (!fs::exists(output_settings[m]))
                    fs::copy_file(input_settings[m], output_settings[m]);
            }
            fs::path previous_seeds = input_dir[MC] / "rng.seed";
            resume[MC] = fs::exists(previous_seeds);
            rng_seed_path = resume[MC] ? previous_seeds : fs::path(pr["s"].as<std::string>());

            resume[MD] = fs::exists(input_velocities);
        }

        std::array<fs::path, 2> input_dir, input_settings{}, output_dir,
                output_settings{output_dir[MC] / "input", output_dir[MD] / "input"},
                input_positions{input_dir[MC] / "positions", input_dir[MD] / "positions"},
                output_positions{output_dir[MC] / "positions", output_dir[MD] / "positions"};

        fs::path input_velocities{input_dir[MD] / "velocities"},
                output_velocities{output_dir[MD] / "velocities"}, rng_seed_path;
        size_t n_bins;
        std::array<bool, 2> sample, resume{};
        bool warmup;
    };
}// namespace ex07
#endif//ESERCIZI_LSN_07_OPTS_HPP
