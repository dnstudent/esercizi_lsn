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

namespace ex09 {

    enum CrossAlgo { Exercise, ExerciseMod, MyAlgo1, MyAlgo2, Fusion, Dummy };

    std::string tag_from(CrossAlgo algo) {
        static std::map<CrossAlgo, std::string> algo_map{{Exercise, "ex"},   {ExerciseMod, "exmod"},
                                                         {MyAlgo1, "my1"},   {MyAlgo2, "my2"},
                                                         {Fusion, "fusion"}, {Dummy, "dummy"}};
        return algo_map[algo];
    }

    class ExOptions {
    private:
        std::map<std::string, CrossAlgo> algo_mapper{{"ex", Exercise},   {"exmod", ExerciseMod},
                                                     {"my1", MyAlgo1},   {"my2", MyAlgo2},
                                                     {"fusion", Fusion}, {"dummy", Dummy}};

    public:
        explicit ExOptions(cxxopts::ParseResult &pr)
            : out_dir(pr["out"].as<fs::path>()), in_path(pr["in"].as<fs::path>()),
              seeds_path(pr["s"].as<std::string>()), primes_path(pr["p"].as<std::string>()),
              primes_line(pr["l"].as<size_t>()), n_iter(pr["n"].as<size_t>()),
              pop_size(pr["m"].as<size_t>()), mut_rate(pr["r"].as<double>()),
              fusion_p(pr["f"].as<double>()), algo(algo_mapper[pr["crossover"].as<std::string>()]) {
            if (!fs::exists(out_dir)) fs::create_directories(out_dir);
            utils::require_existence(in_path);
        }

        const fs::path out_dir, in_path;
        const std::string seeds_path, primes_path;
        const size_t primes_line, n_iter, pop_size;
        const double mut_rate, fusion_p;
        const CrossAlgo algo;
    };
}// namespace ex09

#endif//ESERCIZI_LSN_10_OPTIONS_HPP
