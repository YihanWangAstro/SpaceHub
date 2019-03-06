//
// Created by yihan on 1/3/19.
//

#ifndef SPACEHUB_EOM_H
#define SPACEHUB_EOM_H

#include "type-class.h"
#include "core-computation.h"

namespace SpaceH{

    template <typename Derived, bool IsVelDep>
    class Interactions {
    public:
        SPACEHUB_USING_TYPE_SYSTEM_OF(Derived);

        static constexpr bool is_vel_dep{false};

        template<typename Particles>
        void eval_acc(Particles const &partc, Coord& acc) {
            static_cast<Derived*>(this)->impl_eval_acc(partc, acc);
        }
    private:
        Interactions() = default;
        friend Derived;
    };

    template <typename Derived>
    class Interactions<Derived, true> {
    public:
        SPACEHUB_USING_TYPE_SYSTEM_OF(Derived);

        static constexpr bool is_vel_dep{true};

        template<typename Particles>
        void eval_acc(Particles const &partc, Coord& acc) {
            static_cast<Derived*>(this)->impl_eval_acc(partc, acc);
        }

        template<typename Particles>
        void eval_vel_indep_acc(Particles const &partc, Coord& acc) {
            static_cast<Derived*>(this)->impl_eval_vel_indep_acc(partc, acc);
        }

        template<typename Particles>
        void eval_vel_dep_acc(Particles const &partc, Coord& acc) {
            static_cast<Derived*>(this)->impl_eval_vel_dep_acc(partc, acc);
        }
    private:
        Interactions() = default;
        friend Derived;
    };
}

#endif //SPACEHUB_EOM_H