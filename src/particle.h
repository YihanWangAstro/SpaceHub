
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
        size_t idn{0};
    };

    template <typename Derived>
    class SoAParticles {
    public:
        DECLARE_CRTP_ACCESSOR(mass, auto, Derived);
        DECLARE_CRTP_ACCESSOR(idn, auto, Derived);
        DECLARE_CRTP_ACCESSOR(time, auto, Derived);
        DECLARE_CRTP_ACCESSOR(pos, auto, Derived);
        DECLARE_CRTP_ACCESSOR(vel, auto, Derived);

        void resize(size_t new_sz) {
            static_cast<Derived*>(this)->impl_resize(new_sz);
        }

        void reserve(size_t new_cap) {
            static_cast<Derived*>(this)->impl_reserve(new_cap);
        }

        size_t number() {
            return static_cast<Derived*>(this)->impl_number();
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



    template<typename TypeSystem>
    class SoAPointParticle : public SoAParticles<SoAPointParticle<TypeSystem>> {
    public:
        SPACEHUB_USING_TYPE_SYSTEM_OF(TypeSystem);

        SPACEHUB_STD_ACCESSOR(impl_mass, mass_);
        SPACEHUB_STD_ACCESSOR(impl_idn, idn_);
        SPACEHUB_STD_ACCESSOR(impl_time, time_);
        SPACEHUB_STD_ACCESSOR(impl_pos, pos_);
        SPACEHUB_STD_ACCESSOR(impl_vel, vel_);

        SoAPointParticle() = delete;

        template<typename STL>
        SoAPointParticle(STL const &partc, Scalar t) {
            size_t input_num = partc.size();
            this->reserve(input_num);
            for (auto &p : partc) {
                pos_.emplace_back(p.pos);
                vel_.emplace_back(p.vel);
                mass_.emplace_back(p.mass);
                idn_.emplace_back(p.idn);
            }
            time_ = t;
            active_num = input_num;
        }

        void impl_resize(size_t new_sz) {
            SpaceH::resize_all(new_sz, pos_, vel_, mass_, idn_);
            active_num = new_sz;
        }

        void impl_reserve(size_t new_cap) {
            SpaceH::reserve_all(new_cap, pos_, vel_, mass_, idn_);
        }

        size_t impl_number() {
            return active_num;
        }

    private:
        Coord pos_;
        Coord vel_;
        ScalarArray mass_;
        IdxArray idn_;
        Scalar time_;
        size_t active_num{0};
    };
}
#endif

