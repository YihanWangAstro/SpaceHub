//
// Created by yihan on 1/3/19.
//

#ifndef SPACEHUB_EOM_H
#define SPACEHUB_EOM_H

#include "type_class.h"
#include "core_computation.h"

namespace SpaceH{

    template <typename Derived, bool IsVelDep>
    class Interactions {
    public:
        SPACEHUB_USING_TYPE_SYSTEM_OF(Derived);
        static constexpr bool isVelDependent{false};

        Interactions() {
            calcu::array_set_zeros(ax);
            calcu::array_set_zeros(ay);
            calcu::array_set_zeros(az);
        }

        void resize(size_t new_sz) {
            if constexpr (array_size == SpaceH::DYNAMICAL) {
                SpaceH::resize_all(new_sz, ax, ay, az);
            } else {
                SPACEHUB_ABORT("Fixed size arrays are not allowed to resize!");
            }
        }

        void reserve(size_t new_cap) {
            if constexpr (array_size == SpaceH::DYNAMICAL) {
                SpaceH::reserve_all(new_cap, ax, ay, az);
            } else {
                SPACEHUB_ABORT("Fixed size arrays are not allowed to reserve!");
            }
        }

        template<typename Particles>
        void eval_acc(Particles const &partc) {
            static_cast<Derived*>(this)->impl_eval_acc(partc);
        }

        ScalarArray ax;
        ScalarArray ay;
        ScalarArray az;
    };

    template <typename Derived>
    class Interactions<Derived, true> {
    public:
        SPACEHUB_USING_TYPE_SYSTEM_OF(Derived);
        static constexpr bool isVelDependent{true};

        Interactions() : vel_indep_ax(ax), vel_indep_ay(ay), vel_indep_az(az) {
            calcu::array_set_zeros(ax);
            calcu::array_set_zeros(ay);
            calcu::array_set_zeros(az);
            calcu::array_set_zeros(auxi_vx);
            calcu::array_set_zeros(auxi_vy);
            calcu::array_set_zeros(auxi_vz);
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
            calcu::array_add(ax, vel_indep_ax, vel_dep_ax);
            calcu::array_add(ay, vel_indep_ay, vel_dep_ay);
            calcu::array_add(az, vel_indep_az, vel_dep_az);
        }

        ScalarArray ax;
        ScalarArray ay;
        ScalarArray az;
        ScalarArray vel_dep_ax;
        ScalarArray vel_dep_ay;
        ScalarArray vel_dep_az;
        ScalarArray const &vel_indep_ax;
        ScalarArray const &vel_indep_ay;
        ScalarArray const &vel_indep_az;
        ScalarArray auxi_vx;
        ScalarArray auxi_vy;
        ScalarArray auxi_vz;
    };



    template <typename TypeSystem>
    class NewtonianGrav {
    public:
        SPACEHUB_USING_TYPE_SYSTEM_OF(TypeSystem);
        static constexpr bool isVelDependent{false};

        NewtonianGrav(){
            calcu::array_set_zeros(ax);
            calcu::array_set_zeros(ay);
            calcu::array_set_zeros(az);
        }

        template<typename Particles>
        void eval_acc(Particles const &partc) {
            size_t num = partc.number();
            for (size_t i = 0; i < num; ++i) {
                for (size_t j = i + 1; j < num; ++j) {
                    auto dx  = partc.px[j] - partc.px[i];
                    auto dy  = partc.py[j] - partc.py[i];
                    auto dz  = partc.pz[j] - partc.pz[i];
                    auto r   = sqrt(dx * dx + dy * dy + dz * dz);
                    auto rr3 = 1 / (r * r * r);
                    ax[i] += dx * rr3 * partc.mass[j];
                    ay[i] += dy * rr3 * partc.mass[j];
                    az[i] += dz * rr3 * partc.mass[j];
                    ax[j] -= dx * rr3 * partc.mass[i];
                    ay[j] -= dy * rr3 * partc.mass[i];
                    az[j] -= dz * rr3 * partc.mass[i];
                }
            }
        }

        ScalarArray ax;
        ScalarArray ay;
        ScalarArray az;
    };

}

#endif //SPACEHUB_EOM_H