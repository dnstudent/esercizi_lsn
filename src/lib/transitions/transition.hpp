//
// Created by Davide Nicoli on 25/08/22.
//

#ifndef ESERCIZI_LSN_TRANSITION_HPP
#define ESERCIZI_LSN_TRANSITION_HPP

namespace transitions {
    /*template<typename prob_space, typename state_space, class URBG>
    class IConditionalTransitionOld {
    public:
        virtual ~IConditionalTransitionOld() = default;

        typedef state_space StateSpace;
        typedef prob_space ProbSpace;

        virtual StateSpace sample(URBG &rng, const StateSpace &from) = 0;
        // virtual ProbSpace p(const StateSpace &to, const StateSpace &from) const = 0;
        virtual ProbSpace logp(const StateSpace &to, const StateSpace &from) const = 0;
    };*/

    //    template<class Implementation, typename prob_space, typename state_space>
    //    class IConditionalTransition {
    //    public:
    //        virtual ~IConditionalTransition() = default;
    //
    //        typedef typename Implementation::state_space StateSpace;
    //        typedef typename Implementation::prob_space ProbSpace;
    //
    //        template<class URBG>
    //        StateSpace sample(const StateSpace &from, URBG &rng) {
    //            return static_cast<Implementation *>(this)->_sample(from, rng);
    //        }
    //        ProbSpace logp(const StateSpace &to, const StateSpace &from) const {
    //            return static_cast<const Implementation *>(this)->_logp(to, from);
    //        }
    //    };
}// namespace transitions

#endif//ESERCIZI_LSN_TRANSITION_HPP
