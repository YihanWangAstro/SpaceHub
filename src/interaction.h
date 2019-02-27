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

        static constexpr bool isVelDependent{false};

        template<typename Particles>
        void eval_acc(Particles const &partc, ScalarArray& ax, ScalarArray& ay, ScalarArray& az) {
            static_cast<Derived*>(this)->impl_eval_acc(partc, ax, ay, az);
        }
    };

    template <typename Derived>
    class Interactions<Derived, true> {
    public:
        SPACEHUB_USING_TYPE_SYSTEM_OF(Derived);


        static constexpr bool isVelDependent{true};

        template<typename Particles>
        void eval_acc(Particles const &partc, ScalarArray& ax, ScalarArray& ay, ScalarArray& az) {
            static_cast<Derived*>(this)->impl_eval_acc(partc, ax, ay, az);
        }

        template<typename Particles>
        void eval_vel_indep_acc(Particles const &partc, ScalarArray& ax, ScalarArray& ay, ScalarArray& az) {
            static_cast<Derived*>(this)->impl_eval_vel_indep_acc(partc, ax, ay, az);
        }

        template<typename Particles>
        void eval_vel_dep_acc(Particles const &partc, ScalarArray& ax, ScalarArray& ay, ScalarArray& az) {
            static_cast<Derived*>(this)->impl_eval_vel_dep_acc(partc, ax, ay, az);
        }

        template<typename Particles>
        void eval_aux_vel_dep_acc(Particles const &partc, ScalarArray &ax, ScalarArray &ay, ScalarArray &az) {
            static_cast<Derived*>(this)->impl_eval_aux_vel_dep_acc(partc, ax, ay, az);
        }
    };
}

#endif //SPACEHUB_EOM_H