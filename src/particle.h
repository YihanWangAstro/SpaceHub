
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

        PointParticle(Scalar m, Vector const &p, Vector const &v)
                : pos(p), vel(v), mass(m) {}

        PointParticle(Scalar m, Scalar px, Scalar py, Scalar pz, Scalar vx, Scalar vy, Scalar vz)
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

        friend std::ostream &operator<<(std::ostream &os, SoAParticles const &ps) {
            size_t num = ps.number();
            os << ps.time() << ' ';
            for (size_t i = 0; i < num; ++i) {
                SpaceH::display(os, ps.idn()[i], ps.mass()[i], ps.pos().x[i], ps.pos().y[i], ps.pos().z[i], ps.vel().x[i], ps.vel().y[i], ps.vel().z[i]);
            }
            return os;
        }

    private:
        SoAParticles() = default;
        friend Derived;
    };

#define SPACEHUB_PARTICLE_TYPE_CHECK(CTR, VAL) static_assert(std::is_base_of_v<typename CTR::value_type, VAL>, "Class can only be initialized by containers with its internal 'Particle' type!");
}
#endif

