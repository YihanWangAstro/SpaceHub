//
// Created by yihan on 2/24/19.
//

#ifndef SPACEHUB_PARTICLE_SYSTEM_H
#define SPACEHUB_PARTICLE_SYSTEM_H


#include "dev-tools.h"

namespace SpaceH {

    template<typename Derived>
    class ParticleSystem {
    public:
        SPACEHUB_USING_TYPE_SYSTEM_OF(Derived);

        DECLARE_CRTP_ACCESSOR(mass, ScalarArray, Derived);

        DECLARE_CRTP_ACCESSOR(idn, IdxArray, Derived);

        DECLARE_CRTP_ACCESSOR(time, Scalar, Derived);

        DECLARE_CRTP_ACCESSOR(pos, Coord, Derived);

        DECLARE_CRTP_ACCESSOR(vel, Coord, Derived);

        size_t number() {
            return static_cast<Derived *>(this)->impl_number();
        }

        void advance_time(Scalar dt) {
            static_cast<Derived *>(this)->impl_advance_time(dt);
        }

        void advance_pos(Scalar stepSize, Coord const &velocity) {
            static_cast<Derived *>(this)->impl_advance_pos(stepSize, velocity);
        }

        void advance_vel(Scalar stepSize, Coord const &acceleration) {
            static_cast<Derived *>(this)->impl_advance_vel(stepSize, acceleration);
        }

        void evaluate_acc(Coord &acceleration) {
            static_cast<Derived *>(this)->impl_evaluate_acc(acceleration);
        }

        void drift(Scalar stepSize) {
            static_cast<Derived *>(this)->impl_drift(stepSize);
        }

        void kick(Scalar stepSize) {
            static_cast<Derived *>(this)->impl_kick(stepSize);
        }

        void pre_iter_process() {
            static_cast<Derived *>(this)->impl_pre_iter_process();
        }

        void post_iter_process() {
            static_cast<Derived *>(this)->impl_post_iter_process();
        }

        void impl_pre_iter_process() {}//default implementation

        void impl_post_iter_process() {}//default implementation
    private:
        ParticleSystem() = default;

        friend Derived;
    };
}

#endif //SPACEHUB_PARTICLE_SYSTEM_H
