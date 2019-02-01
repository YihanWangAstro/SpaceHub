

#ifndef SPACEHUB_PARTICLE_H
#define SPACEHUB_PARTICLE_H

#include "dev_tools.h"
#include "vector/vector3.h"

namespace SpaceH {

    template<typename Real>
    struct PointParticle {
    public:
        using Scalar = Real;
        using Vector = Vec3<Scalar>;

        PointParticle() = default;

        PointParticle(size_t id, Scalar m, Vector const &p, Vector const &v)
            : pos(p), vel(v), mass(m), idn(id) {}

        PointParticle(size_t id, Scalar m, Scalar px, Scalar py, Scalar pz, Scalar vx, Scalar vy, Scalar vz)
            : pos(px,py,pz), vel(vx,vy,vz), mass(m), idn(id) {}

        auto basic_info(){
            return std::tie(idn, mass, pos.x, pos.y, pos.z, vel.x, vel.y, vel.z);
        }

        friend std::ostream &operator<<(std::ostream &os, PointParticle const &particle) {
            SpaceH::display(os, particle.idn, particle.mass, particle.pos, particle.vel);
            return os;
        }

        friend std::istream &operator>>(std::istream &is, PointParticle &particle) {
            SpaceH::input(is, particle.idn, particle.mass, particle.pos, particle.vel);
            return is;
        }

        Vector pos;
        Vector vel;
        Scalar mass;
        size_t idn;
    };

    template<typename Real>
    struct InfniteSizeParticle {
    public:
        using Scalar = Real;
        using Vector = Vec3<Scalar>;

        InfniteSizeParticle() = default;

        InfniteSizeParticle(size_t id, Scalar m, Scalar r, Vector const &p, Vector const &v)
                : pos(p), vel(v), mass(m), radius(r), idn(id) {}

        InfniteSizeParticle(size_t id, Scalar m, Scalar r, Scalar px, Scalar py, Scalar pz, Scalar vx, Scalar vy, Scalar vz)
                : pos(px,py,pz), vel(vx,vy,vz), mass(m), radius(r), idn(id) {}

        auto basic_info(){
            return std::tie(idn, mass, pos.x, pos.y, pos.z, vel.x, vel.y, vel.z);
        }

        friend std::ostream &operator<<(std::ostream &os, InfniteSizeParticle const &particle) {
            SpaceH::display(os, particle.idn, particle.mass, particle.radius, particle.pos, particle.vel);
            return os;
        }

        friend std::istream &operator>>(std::istream &is, InfniteSizeParticle &particle) {
            SpaceH::input(is, particle.idn, particle.mass, particle.radius, particle.pos, particle.vel);
            return is;
        }

        Vector pos;
        Vector vel;
        Scalar mass;
        Scalar radius;
        size_t idn;
    };

}
#endif //SPACEHUB_PARTICLE_H
