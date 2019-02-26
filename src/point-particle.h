
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
        size_t idn{};
    };

    template<typename TypeSystem>
    class SoAPointParticle {
    public:
        SPACEHUB_USING_TYPE_SYSTEM_OF(TypeSystem);

        SPACEHUB_STD_ARRAY_INTERFACES(px, px_);
        SPACEHUB_STD_ARRAY_INTERFACES(py, py_);
        SPACEHUB_STD_ARRAY_INTERFACES(pz, pz_);
        SPACEHUB_STD_ARRAY_INTERFACES(vx, vx_);
        SPACEHUB_STD_ARRAY_INTERFACES(vy, vy_);
        SPACEHUB_STD_ARRAY_INTERFACES(vz, vz_);
        SPACEHUB_STD_ARRAY_INTERFACES(mass, mass_);
        SPACEHUB_STD_ARRAY_INTERFACES(idn, idn_);
        SPACEHUB_STD_SCALAR_INTERFACES(time, time_);

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
                std::tie(idn_[i], mass_[i], px_[i], py_[i], pz_[i], vx_[i], vy_[i], vz_[i]) = p.forward_as_tuple();
                i++;
            }

            active_num = input_num;
            time_ = t;
        }

        void resize(size_t new_sz) {
            TYPE_SYSTEM_RESIZE(new_sz, px_, py_, pz_, vx_, vy_, vz_, mass_, idn_);
            active_num = new_sz;
        }

        void reserve(size_t new_cap) {
            TYPE_SYSTEM_RESERVE(new_cap, px_, py_, pz_, vx_, vy_, vz_, mass_, idn_);
        }

        size_t number() {
            return active_num;
        }

        friend std::ostream &operator<<(std::ostream &os, SoAPointParticle const &ps) {
            size_t num = ps.number();
            for (size_t i = 0; i < num; ++i) {
                SpaceH::display(os, ps.idn_[i], ps.mass_[i], ps.px_[i], ps.py_[i], ps.pz_[i], ps.vx_[i], ps.vy_[i], ps.vz_[i]);
            }
            return os;
        }

        friend std::istream &operator>>(std::istream &is, SoAPointParticle &ps) {
            size_t num = ps.number();
            for (size_t i = 0; i < num; ++i) {
                SpaceH::input(is, ps.idn_[i], ps.mass_[i], ps.px_[i], ps.py_[i], ps.pz_[i], ps.vx_[i], ps.vy_[i], ps.vz_[i]);
            }
            return is;
        }

    private:
        ScalarArray px_;
        ScalarArray py_;
        ScalarArray pz_;
        ScalarArray vx_;
        ScalarArray vy_;
        ScalarArray vz_;
        ScalarArray mass_;
        IndexArray idn_;
        Scalar time_;
        size_t active_num{0};
    };

}
#endif

