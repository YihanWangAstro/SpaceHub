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
        DECLARE_CRTP_ACCESSOR(Derived, auto, mass);

        DECLARE_CRTP_ACCESSOR(Derived, auto, idn);

        DECLARE_CRTP_ACCESSOR(Derived, auto, time);

        DECLARE_CRTP_ACCESSOR(Derived, auto, pos);

        DECLARE_CRTP_ACCESSOR(Derived, auto, vel);

        size_t number() const {
            return static_cast<Derived const*>(this)->impl_number();
        }

        template<typename Scalar>
        void advance_time(Scalar dt) {
            static_cast<Derived *>(this)->impl_advance_time(dt);
        }

        template<typename Coord, typename Scalar>
        void advance_pos(Scalar stepSize, Coord const &velocity) {
            static_cast<Derived *>(this)->impl_advance_pos(stepSize, velocity);
        }

        template<typename Coord, typename Scalar>
        void advance_vel(Scalar stepSize, Coord const &acceleration) {
            static_cast<Derived *>(this)->impl_advance_vel(stepSize, acceleration);
        }

        template<typename Coord>
        void evaluate_acc(Coord &acceleration) const {
            static_cast<Derived const*>(this)->impl_evaluate_acc(acceleration);
        }

        template<typename Scalar>
        void drift(Scalar stepSize) {
            static_cast<Derived *>(this)->impl_drift(stepSize);
        }

        template<typename Scalar>
        void kick(Scalar stepSize) {
            static_cast<Derived *>(this)->impl_kick(stepSize);
        }

        void pre_iter_process() {
            static_cast<Derived *>(this)->impl_pre_iter_process();
        }

        void post_iter_process() {
            static_cast<Derived *>(this)->impl_post_iter_process();
        }

        template<typename STL>
        void to_linear_container(STL & stl) {
            static_assert(is_container<STL>::value, "Only STL-like container can be used");
            static_cast<Derived *>(this)->impl_to_linear_container(stl);
        }

        template<typename STL>
        void load_from_linear_container(STL & stl) {
            static_assert(is_container<STL>::value, "Only STL-like container can be used");
            static_cast<Derived *>(this)->impl_load_from_linear_container(stl);
        }
    private:
        ParticleSystem() = default;

        void impl_pre_iter_process() {}//default implementation

        void impl_post_iter_process() {}//default implementation

        friend Derived;
    };
}

#endif //SPACEHUB_PARTICLE_SYSTEM_H
