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
        DECLARE_STD_ARRAY_INTERFACES(ax, Derived);
        DECLARE_STD_ARRAY_INTERFACES(ay, Derived);
        DECLARE_STD_ARRAY_INTERFACES(az, Derived);

        static constexpr bool isVelDependent{false};

        void resize(size_t new_sz) {
            static_cast<Derived*>(this)->impl_resize(new_sz);
        }

        void reserve(size_t new_cap) {
            static_cast<Derived*>(this)->impl_reserve(new_cap);
        }

        template<typename Particles>
        void eval_acc(Particles const &partc) {
            static_cast<Derived*>(this)->impl_eval_acc(partc);
        }
    };

    template <typename Derived>
    class Interactions<Derived, true> {
    public:
        SPACEHUB_USING_TYPE_SYSTEM_OF(Derived);
        DECLARE_STD_ARRAY_INTERFACES(ax, Derived);
        DECLARE_STD_ARRAY_INTERFACES(ay, Derived);
        DECLARE_STD_ARRAY_INTERFACES(az, Derived);
        DECLARE_STD_ARRAY_INTERFACES(vel_dep_ax, Derived);
        DECLARE_STD_ARRAY_INTERFACES(vel_dep_ay, Derived);
        DECLARE_STD_ARRAY_INTERFACES(vel_dep_az, Derived);
        DECLARE_STD_ARRAY_INTERFACES(vel_indep_ax, Derived);
        DECLARE_STD_ARRAY_INTERFACES(vel_indep_ay, Derived);
        DECLARE_STD_ARRAY_INTERFACES(vel_indep_az, Derived);

        static constexpr bool isVelDependent{true};

        void resize(size_t new_sz) {
            static_cast<Derived*>(this)->impl_resize(new_sz);
        }

        void reserve(size_t new_cap) {
            static_cast<Derived*>(this)->impl_reserve(new_cap);
        }

        template<typename Particles>
        void eval_acc(Particles const &partc) {
            static_cast<Derived*>(this)->impl_eval_acc(partc);
        }

        template<typename Particles>
        void eval_vel_indep_acc(Particles const &partc) {
            static_cast<Derived*>(this)->impl_eval_vel_indep_acc(partc);
        }

        template<typename Particles>
        void eval_vel_dep_acc(Particles const &partc) {
            static_cast<Derived*>(this)->impl_eval_vel_dep_acc(partc);
        }

        template<typename Particles>
        void eval_auxi_vel_dep_acc(Particles const &partc) {
            static_cast<Derived*>(this)->impl_eval_auxi_vel_dep_acc(partc);
        }

        void sum_tot_acc() {
            static_cast<Derived*>(this)->impl_sum_tot_acc();
        }
    };
}

#endif //SPACEHUB_EOM_H