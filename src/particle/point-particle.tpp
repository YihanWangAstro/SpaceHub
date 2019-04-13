//
// Created by yihan on 4/3/19.
//

#ifndef SPACEHUB_POINT_PARTICLE_H
#define SPACEHUB_POINT_PARTICLE_H

#include "../particle.h"

namespace space{
    template<typename TypeSystem>
    class SoAPointParticles : public SoAParticles<SoAPointParticles<TypeSystem>> {
    public:
        SPACEHUB_USING_TYPE_SYSTEM_OF(TypeSystem);

        using Particle = PointParticle<Scalar>;

        SPACEHUB_STD_ACCESSOR(auto, impl_mass, mass_);
        SPACEHUB_STD_ACCESSOR(auto, impl_idn, idn_);
        SPACEHUB_STD_ACCESSOR(auto, impl_time, time_);
        SPACEHUB_STD_ACCESSOR(auto, impl_pos, pos_);
        SPACEHUB_STD_ACCESSOR(auto, impl_vel, vel_);

        SoAPointParticles() = delete;

        SoAPointParticles(SoAPointParticles const &)= default;

        SoAPointParticles(SoAPointParticles &&)= default;

        SoAPointParticles &operator=(SoAPointParticles const &)= default;

        SoAPointParticles &operator=(SoAPointParticles &&)= default;

        template<typename STL>
        SoAPointParticles( Scalar t, STL const &partc) {
            static_assert(is_container_v<STL>, "Only STL-like container can be used");
            SPACEHUB_PARTICLE_TYPE_CHECK(STL, Particle);

            size_t input_num = partc.size();
            this->reserve(input_num);
            size_t id = 0;
            for (auto &p : partc) {
                pos_.emplace_back(p.pos);
                vel_.emplace_back(p.vel);
                mass_.emplace_back(p.mass);
                idn_.emplace_back(id++);
            }
            time_ = t;
            active_num = input_num;
        }

        void impl_resize(size_t new_sz) {
            space::resize_all(new_sz, pos_, vel_, mass_, idn_);
            active_num = new_sz;
        }

        void impl_reserve(size_t new_cap) {
            space::reserve_all(new_cap, pos_, vel_, mass_, idn_);
        }

        size_t impl_number() const {
            return active_num;
        }

        friend std::ostream &operator<<(std::ostream &os, SoAPointParticles const &ps) {
            size_t num = ps.number();
            os << ps.time() << ' ';
            for (size_t i = 0; i < num; ++i) {
                space::display(os, ps.idn()[i], ps.mass()[i], ps.pos().x[i], ps.pos().y[i], ps.pos().z[i], ps.vel().x[i], ps.vel().y[i], ps.vel().z[i]);
            }
            return os;
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


#endif //SPACEHUB_POINT_PARTICLE_H
