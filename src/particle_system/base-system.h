
#ifndef GENPARTICLESYSTEM_H
#define GENPARTICLESYSTEM_H

#include "../core_computation.h"
#include "../dev_tools.h"
#include "../macros.h"
#include "../type_class.h"
#include "../particle-system.h"

namespace SpaceH {

    template<typename Particles, typename Interaction>
    class BaseSystem : ParticleSystem<BaseSystem<Particles, Interaction>>{
    public:
        SPACEHUB_USING_TYPE_SYSTEM_OF(Particles);

        BaseSystem() = delete;

        template<typename Container>
        BaseSystem(Container const &partc, Scalar t) : partc_(partc, t) {
            if constexpr (Interaction::isVelDependent) {
                eom_.auxi_vx = partc_.vx;
                eom_.auxi_vy = partc_.vy;
                eom_.auxi_vz = partc_.vz;
            }
        }

        size_t impl_number() {
            return partc_.number();
        }

        void impl_advance_time(Scalar dt) {
            partc_.time_ += dt;
        }

        void impl_advance_pos(Scalar stepSize) {
            advance_pos(partc_.vx, partc_.vy, partc_.vz, stepSize);
        }

        void impl_advance_pos(ScalarArray const &vx, ScalarArray const &vy, ScalarArray const &vz, Scalar stepSize) {
            calcu::array_advance(partc_.px, partc_.vx, stepSize);
            calcu::array_advance(partc_.py, partc_.vy, stepSize);
            calcu::array_advance(partc_.pz, partc_.vz, stepSize);
        }

        void impl_advance_vel(Scalar stepSize) {
            advance_vel(eom_.ax, eom_.ay, eom_.az, stepSize);
        }

        void impl_advance_vel(ScalarArray const &ax, ScalarArray const &ay, ScalarArray const &az, Scalar stepSize) {
            calcu::array_advance(partc_.vx, ax, stepSize);
            calcu::array_advance(partc_.vy, ay, stepSize);
            calcu::array_advance(partc_.vz, az, stepSize);
        }

        void impl_evaluate_acc() {
            eom_.eval_acc(partc_);
        }

        void impl_drift(Scalar stepSize) {
            advance_time(stepSize);
            advance_pos(stepSize);
        }

        void impl_kick(Scalar stepSize) {
            if constexpr (Interaction::isVelDependent) {
                eom_.eval_vel_indep_acc(partc_);

                eom_.eval_vel_dep_acc(partc_);
                eom_.sum_tot_acc();
                calcu::array_advance(eom_.auxi_vx, eom_.ax, 0.5 * stepSize);
                calcu::array_advance(eom_.auxi_vy, eom_.ay, 0.5 * stepSize);
                calcu::array_advance(eom_.auxi_vz, eom_.az, 0.5 * stepSize);

                eom_.eval_auxi_vel_dep_acc(partc_);
                eom_.sum_tot_acc();
                advance_vel(stepSize);

                eom_.eval_vel_dep_acc(partc_);
                eom_.sum_tot_acc();
                calcu::array_advance(eom_.auxi_vx, eom_.ax, 0.5 * stepSize);
                calcu::array_advance(eom_.auxi_vy, eom_.ay, 0.5 * stepSize);
                calcu::array_advance(eom_.auxi_vz, eom_.az, 0.5 * stepSize);
            } else {
                eom_.eval_acc(partc_);
                advance_vel(stepSize);
            }
        }

        friend std::ostream&operator<<(std::ostream& os, BaseSystem const & ps) {
            os << ps.partc_;
        }
    private:
        Particles partc_;
        Interaction eom_;
    };



}

#endif
