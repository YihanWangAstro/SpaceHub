//
// Created by yihan on 1/3/19.
//

#ifndef SPACEHUB_EOM_H
#define SPACEHUB_EOM_H

#include "type_class.h"
#include "core_computation.h"

namespace SpaceH{

    template <typename TypeSystem>
    class NewtonianGrav {
    public:
        SPACEHUB_USING_TYPE_SYSTEM_OF(TypeSystem);
        static constexpr bool isVelDependent{false};

        NewtonianGrav(){
            comp::array_set_zeros(ax);
            comp::array_set_zeros(ay);
            comp::array_set_zeros(az);
        }

        template<typename Particles>
        void eval_acc(const Particles &partc) {
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


    template <typename TypeSystem>
    class PostNewtonianGrav {
    public:
        SPACEHUB_USING_TYPE_SYSTEM_OF(TypeSystem);
        static constexpr bool isVelDependent{true};

        PostNewtonianGrav() : vel_indep_ax(ax), vel_indep_ay(ay), vel_indep_az(az){
            comp::array_set_zeros(ax);
            comp::array_set_zeros(ay);
            comp::array_set_zeros(az);
            comp::array_set_zeros(auxi_vx);
            comp::array_set_zeros(auxi_vy);
            comp::array_set_zeros(auxi_vz);
        }

        template<typename Particles>
        void eval_acc(const Particles &partc) {
            eval_vel_indep_acc(partc);
            eval_vel_dep_acc(partc);
            sum_tot_acc();
        }

        template<typename Particles>
        void eval_vel_indep_acc(const Particles &partc) {
            size_t num = partc.number();
            for (size_t i = 0; i < num; ++i) {
                for (size_t j = i + 1; j < num; ++j) {
                    auto dx  = partc.px[j] - partc.px[i];
                    auto dy  = partc.py[j] - partc.py[i];
                    auto dz  = partc.pz[j] - partc.pz[i];
                    auto r   = sqrt(dx * dx + dy * dy + dz * dz);
                    auto rr3 = 1 / (r * r * r);
                    vel_indep_ax[i] += dx * rr3 * partc.mass[j];
                    vel_indep_ay[i] += dy * rr3 * partc.mass[j];
                    vel_indep_az[i] += dz * rr3 * partc.mass[j];
                    vel_indep_ax[j] -= dx * rr3 * partc.mass[i];
                    vel_indep_ay[j] -= dy * rr3 * partc.mass[i];
                    vel_indep_az[j] -= dz * rr3 * partc.mass[i];
                }
            }
        }

        template<typename Particles>
        void eval_vel_dep_acc(const Particles &partc) {

        }

        template<typename Particles>
        void sum_tot_acc() {

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

}

#endif //SPACEHUB_EOM_H