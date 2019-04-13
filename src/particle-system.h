//
// Created by yihan on 2/24/19.
//

#ifndef SPACEHUB_PARTICLE_SYSTEM_H
#define SPACEHUB_PARTICLE_SYSTEM_H

#include "dev-tools.h"

namespace space {

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
        void advance_pos(Scalar step_size, Coord const &velocity) {
            static_cast<Derived *>(this)->impl_advance_pos(step_size, velocity);
        }

        template<typename Coord, typename Scalar>
        void advance_vel(Scalar step_size, Coord const &acceleration) {
            static_cast<Derived *>(this)->impl_advance_vel(step_size, acceleration);
        }

        template<typename Coord>
        void evaluate_acc(Coord &acceleration) const {
            static_cast<Derived const*>(this)->impl_evaluate_acc(acceleration);
        }

        template<typename Scalar>
        void drift(Scalar step_size) {
            static_cast<Derived *>(this)->impl_drift(step_size);
        }

        template<typename Scalar>
        void kick(Scalar step_size) {
            static_cast<Derived *>(this)->impl_kick(step_size);
        }

        void pre_iter_process() {
            static_cast<Derived *>(this)->impl_pre_iter_process();
        }

        void post_iter_process() {
            static_cast<Derived *>(this)->impl_post_iter_process();
        }

        template<typename STL>
        void to_linear_container(STL & stl) {
            static_assert(is_container_v<STL>, "Only STL-like container can be used");
            static_cast<Derived *>(this)->impl_to_linear_container(stl);
        }

        template<typename STL>
        void load_from_linear_container(STL const& stl) {
            static_assert(is_container_v<STL>, "Only STL-like container can be used");
            static_cast<Derived *>(this)->impl_load_from_linear_container(stl);
        }

        friend std::ostream &operator<<(std::ostream &os, ParticleSystem const &ps) {
            os << static_cast<Derived const&>(ps);
            return os;
        }

        friend std::istream &operator>>(std::istream &is, ParticleSystem &ps) {
            is >> static_cast<Derived&>(ps);
            return is;
        }
    private:
        ParticleSystem() = default;

        void impl_pre_iter_process() {}//default implementation

        void impl_post_iter_process() {}//default implementation

        friend Derived;
    };

    template <typename T>
    constexpr bool is_particle_system_v = std::is_base_of_v<ParticleSystem<T>, T>;
}

#endif //SPACEHUB_PARTICLE_SYSTEM_H
