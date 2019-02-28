
#ifndef GENPARTICLESYSTEM_H
#define GENPARTICLESYSTEM_H

#include "core-computation.h"
#include "dev-tools.h"
#include "../macros.h"
#include "type-class.h"
#include "../particle-system.h"
#include "../accelerations.h"
#include <memory>

namespace SpaceH {

    template<typename Particles, typename Interactions>
    class BaseSystem : ParticleSystem<BaseSystem<Particles, Interactions>> {
    public:
        SPACEHUB_USING_TYPE_SYSTEM_OF(Particles);

        SPACEHUB_STD_ARRAY_INTERFACES(px, partc_.px());
        SPACEHUB_STD_ARRAY_INTERFACES(py, partc_.py());
        SPACEHUB_STD_ARRAY_INTERFACES(pz, partc_.pz());
        SPACEHUB_STD_ARRAY_INTERFACES(vx, partc_.vx());
        SPACEHUB_STD_ARRAY_INTERFACES(vy, partc_.vy());
        SPACEHUB_STD_ARRAY_INTERFACES(vz, partc_.vz());
        SPACEHUB_STD_ARRAY_INTERFACES(mass, partc_.mass());
        SPACEHUB_STD_ARRAY_INTERFACES(idn, partc_.idn());
        SPACEHUB_STD_ARRAY_INTERFACES(ax, acc_.ax());
        SPACEHUB_STD_ARRAY_INTERFACES(ay, acc_.ay());
        SPACEHUB_STD_ARRAY_INTERFACES(az, acc_.az());

        BaseSystem() = delete;

        template<typename Container>
        BaseSystem(Container const &partc, Scalar t) : partc_(partc, t), acc_(partc.size()) {}

        size_t impl_number() {
            return partc_.number();
        }

        void impl_advance_time(Scalar dt) {
            partc_.time() += dt;
        }

        void impl_advance_pos(Scalar stepSize) {
            impl_advance_pos(vx(), vy(), vz(), stepSize);
        }

        void impl_advance_pos(ScalarArray const &vx, ScalarArray const &vy, ScalarArray const &vz, Scalar stepSize) {
            calc::array_advance(px(), vx, stepSize);
            calc::array_advance(py(), vy, stepSize);
            calc::array_advance(pz(), vz, stepSize);
        }

        void impl_advance_vel(Scalar stepSize) {
            impl_advance_vel(ax(), ay(), az(), stepSize);
        }

        void impl_advance_vel(ScalarArray const &a_x, ScalarArray const &a_y, ScalarArray const &a_z, Scalar stepSize) {
            calc::array_advance(vx(), a_x, stepSize);
            calc::array_advance(vy(), a_y, stepSize);
            calc::array_advance(vz(), a_z, stepSize);
        }

        void impl_evaluate_acc(ScalarArray &a_x, ScalarArray &a_y, ScalarArray &a_z) {
            eom_.eval_acc(partc_, a_x, a_y, a_z);
        }

        void impl_drift(Scalar stepSize) {
            impl_advance_time(stepSize);
            impl_advance_pos(stepSize);
        }

        void impl_kick(Scalar stepSize) {
            if constexpr (Interactions::isVelDependent) {
                Scalar halfStep = 0.5 * stepSize;
                eom_.eval_vel_indep_acc(partc_, acc_.vid_ax(), acc_.vid_ay(), acc_.vid_az());
                kick_pseu_vel(halfStep);
                kick_real_vel(stepSize);
                kick_pseu_vel(halfStep);
            } else {
                eom_.eval_acc(partc_, ax(), ay(), az());
                impl_advance_vel(stepSize);
            }
        }

        friend std::ostream &operator<<(std::ostream &os, BaseSystem const &ps) {
            os << ps.partc_;
        }

    private:
        void kick_pseu_vel(Scalar stepSize){
            eom_.eval_vel_dep_acc(partc_, acc_.vd_ax(), acc_.vd_ay(), acc_.vd_az());
            sum_all_acc(acc_);
            calc::array_advance(partc_.aux_vx(), ax(), stepSize);
            calc::array_advance(partc_.aux_vy(), ay(), stepSize);
            calc::array_advance(partc_.aux_vz(), az(), stepSize);
        }

        void kick_real_vel(Scalar stepSize){
            eom_.eval_aux_vel_dep_acc(partc_, acc_.vd_ax(), acc_.vd_ay(), acc_.vd_az());
            sum_all_acc(acc_);
            impl_advance_vel(stepSize);
        }

        Particles partc_;
        Interactions eom_;
        Accelerations<ScalarArray, Interactions::isVelDependent> acc_;
    };
}

#endif
