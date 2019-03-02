//
// Created by yihan on 2/25/19.
//

#ifndef SPACEHUB_ARCHAIN_H
#define SPACEHUB_ARCHAIN_H

#include "regu-system.h"

namespace SpaceH {

    template<typename Particles, typename Interactions>
    class ChainSystem : public ParticleSystem<ChainSystem<Particles, Interactions>> {
    public:
        SPACEHUB_USING_TYPE_SYSTEM_OF(Particles);

        SPACEHUB_STD_INTERFACES(mass, ptc_.mass());

        SPACEHUB_STD_INTERFACES(idn, ptc_.idn());

        SPACEHUB_STD_INTERFACES(pos, ptc_.pos());

        SPACEHUB_STD_INTERFACES(vel, ptc_.vel());

        SPACEHUB_STD_INTERFACES(time, ptc_.time());

        SPACEHUB_STD_INTERFACES(acc, acc_.acc());

        ChainSystem() = delete;

        template<typename Container>
        ChainSystem(Container const &partc, Scalar t) : ptc_(partc, t), acc_(partc.size()) {}

        size_t impl_number() {
            return ptc_.number();
        }

        void impl_advance_time(Scalar dt) {
            ptc_.time() += dt;
        }

        void impl_advance_pos(Scalar stepSize) {
            impl_advance_pos(vel(), stepSize);
        }

        void impl_advance_pos(Coord const &velocity, Scalar stepSize) {
            calc::coord_advance(pos(), velocity, stepSize);
        }

        void impl_advance_vel(Scalar stepSize) {
            impl_advance_vel(acc(), stepSize);
        }

        void impl_advance_vel(Coord const &acceleration, Scalar stepSize) {
            calc::coord_advance(vel(), acceleration, stepSize);
        }

        void impl_evaluate_acc(Coord &acceleration) {
            eom_.eval_acc(ptc_, acceleration);
        }

        void impl_drift(Scalar stepSize) {
            impl_advance_time(stepSize);
            impl_advance_pos(stepSize);
        }

        void impl_kick(Scalar stepSize) {
            if constexpr (Interactions::is_vel_dep) {
                Scalar halfStep = 0.5 * stepSize;
                eom_.eval_vel_indep_acc(ptc_, acc_.vel_indep_acc());
                kick_pseu_vel(halfStep);
                kick_real_vel(stepSize);
                kick_pseu_vel(halfStep);
            } else {
                eom_.eval_acc(ptc_, acc());
                impl_advance_vel(stepSize);
            }
        }

        friend std::ostream &operator<<(std::ostream &os, ChainSystem const &ps) {
            os << ps.ptc_;
        }

    private:
        void kick_pseu_vel(Scalar stepSize) {
            eom_.eval_vel_dep_acc(ptc_, acc_.vel_dep_acc());
            calc::array_add(acc(), acc_.vel_indep_acc(), acc_.vel_dep_acc());
            calc::coord_advance(ptc_.aux_vel(), acc(), stepSize);
        }

        void kick_real_vel(Scalar stepSize) {
            eom_.eval_aux_vel_dep_acc(ptc_, acc_.vel_dep_acc());
            calc::array_add(acc(), acc_.vel_indep_acc(), acc_.vel_dep_acc());
            impl_advance_vel(stepSize);
        }

        Particles ptc_;
        Interactions eom_;
        Accelerations<ScalarArray, Interactions::is_vel_dep> acc_;
    };

}
#endif //SPACEHUB_ARCHAIN_H
