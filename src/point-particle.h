
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
    class SoAParticle {
    public:
        SPACEHUB_USING_TYPE_SYSTEM_OF(Derived);

        SPACEHUB_STD_INTERFACES(mass, mass_);
        SPACEHUB_STD_INTERFACES(idn, idn_);
        SPACEHUB_STD_INTERFACES(time, time_);
        SPACEHUB_STD_INTERFACES(pos, pos_);
        SPACEHUB_STD_INTERFACES(vel, vel_);

        SoAParticle() = delete;

        template<typename Container>
        SoAParticle(Container const &partc, Scalar t) {
            size_t input_num = partc.size();
            this->reserve(input_num);
            for (auto &p : partc) {
                pos_.emplace_back(p.pos);
                vel_.emplace_back(p.vel);
                mass_.emplace_back(p.mass);
                idn_.emplace_back(p.idn);
            }
            time_ = t;
        }

        void resize(size_t new_sz) {
            static_cast<Derived*>(this)->impl_resize(new_sz);
        }

        void reserve(size_t new_cap) {
            static_cast<Derived*>(this)->impl_reserve(new_cap);
        }

        size_t number() {
            return static_cast<Derived*>(this)->impl_number();
        }

    private:
        Coord pos_;
        Coord vel_;
        ScalarArray mass_;
        IdxArray idn_;
        Scalar time_;
    };

    template<typename TypeSystem, bool IsVelDep>
    class SoAPointParticle {
    public:
        SPACEHUB_USING_TYPE_SYSTEM_OF(TypeSystem);

        SPACEHUB_STD_INTERFACES(mass, mass_);
        SPACEHUB_STD_INTERFACES(idn, idn_);
        SPACEHUB_STD_INTERFACES(time, time_);
        SPACEHUB_STD_INTERFACES(pos, pos_);
        SPACEHUB_STD_INTERFACES(vel, vel_);

        SoAPointParticle() = delete;

        template<typename Container>
        SoAPointParticle(Container const &partc, Scalar t) {
            size_t input_num = partc.size();
            this->reserve(input_num);
            for (auto &p : partc) {
                pos_.emplace_back(p.pos);
                vel_.emplace_back(p.vel);
                mass_.emplace_back(p.mass);
                idn_.emplace_back(p.idn);
                active_num++;
            }
            time_ = t;
        }

        void resize(size_t new_sz) {
            SpaceH::resize_all(new_sz, pos_, vel_, mass_, idn_);
            active_num = new_sz;
        }

        void reserve(size_t new_cap) {
            SpaceH::reserve_all(new_cap, pos_, vel_, mass_, idn_);
        }

        size_t number() {
            return active_num;
        }

        friend std::ostream &operator<<(std::ostream &os, SoAPointParticle const &ps) {
            size_t num = ps.number();
            os << ps.time_ << ' ';
            for (size_t i = 0; i < num; ++i) {
                SpaceH::display(os, ps.idn_[i], ps.mass_[i], ps.pos_.x[i], ps.pos_.y[i], ps.pos_.z[i], ps.vel_.x[i], ps.vel_.y[i], ps.vel_.z[i]);
            }
            return os;
        }

        friend std::istream &operator>>(std::istream &is, SoAPointParticle &ps) {
            size_t num = ps.number();
            is >> ps.time_ << ' ';
            for (size_t i = 0; i < num; ++i) {
                SpaceH::input(is, ps.idn_[i], ps.mass_[i], ps.pos_.x[i], ps.pos_.y[i], ps.pos_.z[i], ps.vel_.x[i], ps.vel_.y[i], ps.vel_.z[i]);
            }
            return is;
        }

    private:
        Coord pos_;
        Coord vel_;
        ScalarArray mass_;
        IdxArray idn_;
        Scalar time_;
        size_t active_num{0};
    };

    template<typename TypeSystem>
    class SoAPointParticle<TypeSystem, true> {
    public:
        SPACEHUB_USING_TYPE_SYSTEM_OF(TypeSystem);

        SPACEHUB_STD_INTERFACES(mass, mass_);
        SPACEHUB_STD_INTERFACES(idn, idn_);
        SPACEHUB_STD_INTERFACES(time, time_);
        SPACEHUB_STD_INTERFACES(pos, pos_);
        SPACEHUB_STD_INTERFACES(vel, vel_);
        SPACEHUB_STD_INTERFACES(aux_vel, aux_vel_);

        SoAPointParticle() = delete;

        template<typename Container>
        SoAPointParticle(Container const &partc, Scalar t) {
            size_t input_num = partc.size();
            this->reserve(input_num);
            for (auto &p : partc) {
                pos_.emplace_back(p.pos);
                vel_.emplace_back(p.vel);
                mass_.emplace_back(p.mass);
                idn_.emplace_back(p.idn);
                active_num++;
            }
            aux_vel_ = vel_;
            time_ = t;
        }

        void resize(size_t new_sz) {
            SpaceH::resize_all(new_sz, pos_, vel_, aux_vel_, mass_, idn_);
            active_num = new_sz;
        }

        void reserve(size_t new_cap) {
            SpaceH::reserve_all(new_cap, pos_, vel_, aux_vel_, mass_, idn_);
        }

        size_t number() {
            return active_num;
        }

        friend std::ostream &operator<<(std::ostream &os, SoAPointParticle const &ps) {
            size_t num = ps.number();
            os << ps.time_ << ' ';
            for (size_t i = 0; i < num; ++i) {
                SpaceH::display(os, ps.idn_[i], ps.mass_[i], ps.pos_.x[i], ps.pos_.y[i], ps.pos_.z[i], ps.vel_.x[i], ps.vel_.y[i], ps.vel_.z[i]);
            }
            return os;
        }

        friend std::istream &operator>>(std::istream &is, SoAPointParticle &ps) {
            size_t num = ps.number();
            is >> ps.time_ << ' ';
            for (size_t i = 0; i < num; ++i) {
                SpaceH::input(is, ps.idn_[i], ps.mass_[i], ps.pos_.x[i], ps.pos_.y[i], ps.pos_.z[i], ps.vel_.x[i], ps.vel_.y[i], ps.vel_.z[i]);
            }
            return is;
        }

    private:
        Coord pos_;
        Coord vel_;
        Coord aux_vel_;
        ScalarArray mass_;
        IdxArray idn_;
        Scalar time_;
        size_t active_num{0};
    };

}
#endif

