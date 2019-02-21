
#ifndef PARTICLES_H
#define PARTICLES_H

#include "dev_tools.h"
#include "vector/vector3.h"
#include "type_class.h"

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
                : pos(px, py, pz), vel(vx, vy, vz), mass(m), idn(id) {}

        auto forward_as_tuple() {
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

    template<typename TypeSystem>
    class SoAPointParticle {
    public:
        SPACEHUB_USING_TYPE_SYSTEM_OF(TypeSystem);

        SoAPointParticle() = delete;

        template<typename Container>
        SoAPointParticle(Container const &partc, Scalar t) {
            size_t input_num = partc.size();

            if constexpr(array_size == SpaceH::DYNAMICAL) {
                this->resize(input_num);
            } else {
                if (input_num > array_size) {
                    SPACEHUB_ABORT("Input particle number exceeds the capacity!");
                }
            }

            size_t i = 0;
            for (auto &p : partc) {
                std::tie(idn[i], mass[i], px[i], py[i], pz[i], vx[i], vy[i], vz[i]) = p.forward_as_tuple();
                i++;
            }

            active_num = input_num;
            time = t;
        }

        void resize(size_t new_sz) {
            if constexpr (array_size == SpaceH::DYNAMICAL) {
                SpaceH::resize_all(new_sz, px, py, pz, vx, vy, vz, mass, idn);
                active_num = new_sz;
            } else {
                SPACEHUB_ABORT("Fixed size arrays are not allowed to resize!");
            }
        }

        void reserve(size_t new_cap) {
            if constexpr (array_size == SpaceH::DYNAMICAL) {
                SpaceH::reserve_all(new_cap, px, py, pz, vx, vy, vz, mass, idn);
            } else {
                SPACEHUB_ABORT("Fixed size arrays are not allowed to reserve!");
            }
        }

        size_t number() {
            return active_num;
        }

        friend std::ostream &operator<<(std::ostream &os, SoAPointParticle const &ps) {
            size_t num = ps.number();
            for (size_t i = 0; i < num; ++i) {
                SpaceH::display(os, ps.idn[i], ps.mass[i], ps.px[i], ps.py[i], ps.pz[i], ps.vx[i], ps.vy[i], ps.vz[i]);
            }
            return os;
        }

        friend std::istream &operator>>(std::istream &is, SoAPointParticle &ps) {
            size_t num = ps.number();
            for (size_t i = 0; i < num; ++i) {
                SpaceH::input(is, ps.idn[i], ps.mass[i], ps.px[i], ps.py[i], ps.pz[i], ps.vx[i], ps.vy[i], ps.vz[i]);
            }
            return is;
        }

        ScalarArray px;
        ScalarArray py;
        ScalarArray pz;
        ScalarArray vx;
        ScalarArray vy;
        ScalarArray vz;
        ScalarArray mass;
        IndexArray idn;
        Scalar time;
    private:
        size_t active_num{0};
    };

}
#endif

