//
// Created by Davide Nicoli on 27/10/22.
//

#ifndef ESERCIZI_LSN_PDF_HPP
#define ESERCIZI_LSN_PDF_HPP

#include <cmath>

namespace distributions {
    template<typename prob_space, typename state_space>
    struct IPDF {
        typedef prob_space ProbSpace;
        typedef state_space StateSpace;

        virtual prob_space logp(const state_space &) const = 0;
        virtual prob_space p(const state_space &x) const { return std::exp(logp(x)); }
    };
}// namespace distributions

#endif//ESERCIZI_LSN_PDF_HPP
