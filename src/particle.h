
#ifndef PARTICLES_H
#define PARTICLES_H

#include "dev-tools.h"
#include "vector/vector3.h"
#include "type-class.h"

namespace SpaceH {

    template<typename Real>
    struct PointParticle {
    public:
        using Scalar = Real;
        using Vector = Vec3<Scalar>;

        PointParticle() = default;

        explicit PointParticle(Scalar m, Vector p, Vector v)
                : pos(p), vel(v), mass(m) {}

        explicit PointParticle(Scalar m, Scalar px = 0, Scalar py = 0, Scalar pz = 0, Scalar vx = 0, Scalar vy = 0, Scalar vz = 0)
                : pos(px, py, pz), vel(vx, vy, vz), mass(m) {}

        friend std::ostream &operator<<(std::ostream &os, PointParticle const &particle) {
            SpaceH::display(os, particle.mass, particle.pos, particle.vel);
            return os;
        }

        friend std::istream &operator>>(std::istream &is, PointParticle &particle) {
            SpaceH::input(is, particle.mass, particle.pos, particle.vel);
            return is;
        }

        Vector pos;
        Vector vel;
        Scalar mass;
    };

    template <typename Derived>
    class SoAParticles {
    public:
        DECLARE_CRTP_ACCESSOR(Derived, auto, mass);
        DECLARE_CRTP_ACCESSOR(Derived, auto, idn);
        DECLARE_CRTP_ACCESSOR(Derived, auto, time);
        DECLARE_CRTP_ACCESSOR(Derived, auto, pos);
        DECLARE_CRTP_ACCESSOR(Derived, auto, vel);

        void resize(size_t new_sz) {
            static_cast<Derived*>(this)->impl_resize(new_sz);
        }

        void reserve(size_t new_cap) {
            static_cast<Derived*>(this)->impl_reserve(new_cap);
        }

        size_t number() const {
            return static_cast<Derived const*>(this)->impl_number();
        }
    private:
        SoAParticles() = default;
        friend Derived;
    };

#define SPACEHUB_PARTICLE_TYPE_CHECK(CTR, VAL) static_assert(std::is_base_of_v<typename CTR::value_type, VAL>, "Class can only be initialized by containers with its internal 'Particle' type!");
}
#endif

