//
// Created by Davide Nicoli on 18/11/22.
//

#ifndef ESERCIZI_LSN_10_OPTIONS_HPP
#define ESERCIZI_LSN_10_OPTIONS_HPP


#include <filesystem>
#include <map>

#include <cxxopts.hpp>

#include "utils.hpp"

namespace fs = std::filesystem;

namespace ex10 {

    enum CrossAlgo { Exercise, ExerciseMod, MyAlgo2, Fusion, Dummy };

    std::string tag_from(CrossAlgo algo) {
        static std::map<CrossAlgo, std::string> algo_map{{Exercise, "ex"},
                                                         {ExerciseMod, "exmod"},
                                                         {MyAlgo2, "my2"},
                                                         {Fusion, "fusion"},
                                                         {Dummy, "dummy"}};
        return algo_map[algo];
    }

    class ExOptions {
    private:
        std::map<std::string, CrossAlgo> algo_mapper{{"ex", Exercise},
                                                     {"exmod", ExerciseMod},
                                                     {"my2", MyAlgo2},
                                                     {"fusion", Fusion},
                                                     {"dummy", Dummy}};

    public:
        explicit ExOptions(cxxopts::ParseResult &pr)
            : out_dir(pr["out"].as<fs::path>()), in_path(pr["in"].as<fs::path>()),
              seeds_path(pr["s"].as<std::string>()), primes_path(pr["p"].as<std::string>()),
              primes_line(pr["l"].as<size_t>()), pop_size(pr["m"].as<size_t>()),
              migration_length(pr["migration_length"].as<size_t>()),
              n_migrations(pr["n_migrations"].as<size_t>()), mut_rate(pr["r"].as<double>()),
              fusion_p(pr["f"].as<double>()), algo(algo_mapper[pr["crossover"].as<std::string>()]) {
            if (!fs::exists(out_dir)) fs::create_directories(out_dir);
            utils::require_existence(in_path);
        }

        const fs::path out_dir, in_path;
        const std::string seeds_path, primes_path;
        const size_t primes_line, pop_size, migration_length, n_migrations;
        const double mut_rate, fusion_p;
        const CrossAlgo algo;
    };
}// namespace ex10

#endif//ESERCIZI_LSN_10_OPTIONS_HPP
