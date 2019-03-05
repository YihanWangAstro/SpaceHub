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

        size_t number() {
            return static_cast<Derived *>(this)->impl_number();
        }

        void advance_time(Scalar dt) {
            static_cast<Derived *>(this)->impl_advance_time(dt);
        }

        void advance_pos(Scalar stepSize) {
            static_cast<Derived *>(this)->impl_advance_pos(stepSize);
        }

        void advance_pos(Coord const &velocity, Scalar stepSize) {
            static_cast<Derived *>(this)->impl_advance_pos(velocity, stepSize);
        }

        void advance_vel(Scalar stepSize) {
            static_cast<Derived *>(this)->impl_advance_vel(stepSize);
        }

        void advance_vel(Coord const &acceleration, Scalar stepSize) {
            static_cast<Derived *>(this)->impl_advance_vel(acceleration, stepSize);
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
    };
}

#endif //SPACEHUB_PARTICLE_SYSTEM_H
