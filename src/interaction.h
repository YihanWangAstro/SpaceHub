//
// Created by yihan on 1/3/19.
//

#ifndef SPACEHUB_EOM_H
#define SPACEHUB_EOM_H

#include "type-class.h"
#include "core-computation.tpp"

namespace space{

    template <typename Derived, bool ExtVelDep, bool ExtVelIndep>
    class Interactions {
    public:
        static constexpr bool has_extra_vel_dep_acc{ExtVelDep};
        static constexpr bool has_extra_vel_indep_acc{ExtVelIndep};

        template<typename Particles>
        void eval_acc(Particles const &partc, typename Particles::Coord& acc) {
            static_cast<Derived*>(this)->impl_eval_acc(partc, acc);
        }

        template<typename Particles>
        void eval_newtonian_acc(Particles const &partc, typename Particles::Coord& acc) {
            static_cast<Derived*>(this)->impl_eval_newtonian_acc(partc, acc);
        }

        template<typename Particles>
        void eval_extra_vel_indep_acc(Particles const &partc, typename Particles::Coord& acc) {
            static_cast<Derived*>(this)->impl_eval_extra_vel_indep_acc(partc, acc);
        }

        template<typename Particles>
        void eval_extra_vel_dep_acc(Particles const &partc, typename Particles::Coord& acc) {
            static_cast<Derived*>(this)->impl_eval_extra_vel_dep_acc(partc, acc);
        }
    private:
        Interactions() = default;
        friend Derived;
    };

    template <typename T, bool a, bool b>
    constexpr bool is_interactions_v = std::is_base_of_v<Interactions<T,a,b>, T>;
}

#endif //SPACEHUB_EOM_H