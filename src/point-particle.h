
#ifndef PARTICLES_H
#define PARTICLES_H

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

        auto forward_as_tuple(){
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

        SPACEHUB_STD_INTERFACES_FOR_ARRAY(px, ScalarArray, px_);

        SPACEHUB_STD_INTERFACES_FOR_ARRAY(py, ScalarArray, py_);

        SPACEHUB_STD_INTERFACES_FOR_ARRAY(pz, ScalarArray, pz_);

        SPACEHUB_STD_INTERFACES_FOR_ARRAY(vx, ScalarArray, vx_);

        SPACEHUB_STD_INTERFACES_FOR_ARRAY(vy, ScalarArray, vy_);

        SPACEHUB_STD_INTERFACES_FOR_ARRAY(vz, ScalarArray, vz_);

        SPACEHUB_STD_INTERFACES_FOR_ARRAY(mass, ScalarArray, mass_);

        SPACEHUB_STD_INTERFACES_FOR_ARRAY(idn, IndexArray, idn_);

        SoAPointParticle() = default;

        template<typename Container>
        explicit SoAPointParticle(Container const &partc) {
            size_t input_num = partc.size();
            if (input_num > _capacity_) {
                SPACEHUB_ABORT("Input particle number exceeds the capacity!");
            } else {
                size_t i = 0;
                for(auto& p : partc){
                    std::tie(idn(i), mass(i), px(i), py(i), pz(i), vx(i), vy(i), vz(i)) = p.forward_as_tuple();
                    i++;
                }
                active_num = input_num;
            }
        }

        size_t number() {
            return active_num;
        }

        friend std::ostream &operator<<(std::ostream &os, SoAPointParticle const &ps) {
            size_t num = ps.number();
            for (size_t i = 0; i < num; ++i) {
                SpaceH::display(os, ps.idn(i), ps.mass(i), ps.px(i), ps.py(i), ps.pz(i), ps.vx(i), ps.vy(i), ps.vz(i));
            }
            return os;
        }

        friend std::istream &operator>>(std::istream &is, SoAPointParticle &ps) {
            size_t num = ps.number();
            for (size_t i = 0; i < num; ++i) {
                SpaceH::input(is, ps.idn(i), ps.mass(i), ps.px(i), ps.py(i), ps.pz(i), ps.vx(i), ps.vy(i), ps.vz(i));
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
        size_t active_num{_capacity_};
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
#endif

