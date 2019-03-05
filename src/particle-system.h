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

        void advance_pos(ScalarArray const &vx, ScalarArray const &vy, ScalarArray const &vz, Scalar stepSize) {
            static_cast<Derived *>(this)->impl_advance_pos(vx, vy, vz, stepSize);
        }

        void advance_vel(Scalar stepSize) {
            static_cast<Derived *>(this)->impl_advance_vel(stepSize);
        }

        void advance_vel(ScalarArray const &ax, ScalarArray const &ay, ScalarArray const &az, Scalar stepSize) {
            static_cast<Derived *>(this)->impl_advance_vel(ax, ay, az, stepSize);
        }

        void evaluate_acc(ScalarArray &ax, ScalarArray &ay, ScalarArray &az) {
            static_cast<Derived *>(this)->impl_evaluate_acc(ax, ay, az);
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
