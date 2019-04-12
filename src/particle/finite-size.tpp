//
// Created by yihan on 4/3/19.
//

#ifndef SPACEHUB_FINITE_SIZE_H
#define SPACEHUB_FINITE_SIZE_H
#include "../particle.h"

namespace SpaceH{
    template<typename TypeSystem>
    class SoAFiniteSizeParticles : public SoAParticles<SoAFiniteSizeParticles<TypeSystem>> {
    public:
        SPACEHUB_USING_TYPE_SYSTEM_OF(TypeSystem);

        struct Particle : public PointParticle<Scalar> {
        public:
            using Scalar = typename TypeSystem::Scalar;
            using Vector = Vec3<Scalar>;

            Particle() = default;

            Particle(Scalar m, Scalar r, Vector p, Vector v)
                    : PointParticle<Scalar>{m, p, v}, radius{r}  {}

            Particle(Scalar m, Scalar r, Scalar px = 0, Scalar py = 0, Scalar pz = 0, Scalar vx = 0, Scalar vy = 0, Scalar vz = 0)
                    : PointParticle<Scalar>{m, px, py, pz, vx, vy, vz}, radius{r} {}

            friend std::ostream &operator<<(std::ostream &os, Particle const &particle) {
                SpaceH::display(os, particle.mass, particle.radius, particle.pos, particle.vel);
                return os;
            }

            friend std::istream &operator>>(std::istream &is, Particle &particle) {
                SpaceH::input(is, particle.mass, particle.radius, particle.pos, particle.vel);
                return is;
            }
            Scalar radius;
        };

        SPACEHUB_STD_ACCESSOR(auto, impl_mass, mass_);

        SPACEHUB_STD_ACCESSOR(auto, impl_idn, idn_);

        SPACEHUB_STD_ACCESSOR(auto, impl_time, time_);

        SPACEHUB_STD_ACCESSOR(auto, impl_pos, pos_);

        SPACEHUB_STD_ACCESSOR(auto, impl_vel, vel_);

        SPACEHUB_STD_ACCESSOR(auto, radius, radius_);

        SoAFiniteSizeParticles() = delete;

        SoAFiniteSizeParticles(SoAFiniteSizeParticles const &) = default;

        SoAFiniteSizeParticles(SoAFiniteSizeParticles &&)= default;

        SoAFiniteSizeParticles &operator=(SoAFiniteSizeParticles const &) = default;

        SoAFiniteSizeParticles &operator=(SoAFiniteSizeParticles &&) = default;

        template<typename STL>
        SoAFiniteSizeParticles(Scalar t, STL const &partc) {
            static_assert(is_container_v<STL>, "Only STL-like container can be used");
            SPACEHUB_PARTICLE_TYPE_CHECK(STL, Particle);

            size_t input_num = partc.size();
            this->reserve(input_num);
            size_t id = 0;
            for (auto &p : partc) {
                pos_.emplace_back(p.pos);
                vel_.emplace_back(p.vel);
                mass_.emplace_back(p.mass);
                radius_.emplace_back(p.radius);
                idn_.emplace_back(id++);
            }
            time_ = t;
            active_num = input_num;
        }

        void impl_resize(size_t new_sz) {
            SpaceH::resize_all(new_sz, pos_, vel_, mass_, radius_, idn_);
            active_num = new_sz;
        }

        void impl_reserve(size_t new_cap) {
            SpaceH::reserve_all(new_cap, pos_, vel_, mass_, radius_, idn_);
        }

        size_t impl_number() const {
            return active_num;
        }

        friend std::ostream &operator<<(std::ostream &os, SoAFiniteSizeParticles const &ps) {
            size_t num = ps.number();
            os << ps.time() << ' ';
            for (size_t i = 0; i < num; ++i) {
                SpaceH::display(os, ps.idn()[i], ps.mass()[i], ps.radius()[i], ps.pos().x[i], ps.pos().y[i], ps.pos().z[i], ps.vel().x[i], ps.vel().y[i], ps.vel().z[i]);
            }
            return os;
        }
    private:
        Coord pos_;
        Coord vel_;
        ScalarArray mass_;
        ScalarArray radius_;
        IdxArray idn_;
        Scalar time_;
        size_t active_num{0};
    };
}
#endif //SPACEHUB_FINITE_SIZE_H
