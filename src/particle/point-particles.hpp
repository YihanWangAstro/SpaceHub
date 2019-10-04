//
// Created by yihan on 4/3/19.
//

#ifndef SPACEHUB_POINT_PARTICLE_H
#define SPACEHUB_POINT_PARTICLE_H

#include "../particle.hpp"

namespace space {
/*---------------------------------------------------------------------------*\
    Class PointParticles Declaration
\*---------------------------------------------------------------------------*/
template <typename TypeSystem>
class PointParticles : public Particles<PointParticles<TypeSystem>> {
 public:
  // Type members
  SPACEHUB_USING_TYPE_SYSTEM_OF(TypeSystem);

  using Base = Particles<PointParticles<TypeSystem>>;

  using Particle = PointParticle<Scalar>;

  // Constructors
  PointParticles() = delete;

  PointParticles(PointParticles const &) = default;

  PointParticles(PointParticles &&) noexcept = default;

  PointParticles &operator=(PointParticles const &) = default;

  PointParticles &operator=(PointParticles &&) noexcept = default;

  template <typename STL>
  PointParticles(Scalar t, STL const &partc);

  CRTP_impl :
      // CRTP implementation

      SPACEHUB_STD_ACCESSOR(auto, impl_mass, mass_);

  SPACEHUB_STD_ACCESSOR(auto, impl_idn, idn_);

  SPACEHUB_STD_ACCESSOR(auto, impl_time, time_);

  SPACEHUB_STD_ACCESSOR(auto, impl_pos, pos_);

  SPACEHUB_STD_ACCESSOR(auto, impl_vel, vel_);

  void impl_resize(size_t new_sz);

  void impl_reserve(size_t new_cap);

  size_t impl_number() const;

 private:
  // Private members
  Coord pos_;

  Coord vel_;

  ScalarArray mass_;

  IdxArray idn_;

  Scalar time_;

  size_t active_num{0};
};

/*---------------------------------------------------------------------------*\
    Class PointParticles Implementation
\*---------------------------------------------------------------------------*/
template <typename TypeSystem>
template <typename STL>
PointParticles<TypeSystem>::PointParticles(Scalar t, const STL &partc) {
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

template <typename TypeSystem>
void PointParticles<TypeSystem>::impl_resize(size_t new_sz) {
  space::resize_all(new_sz, pos_, vel_, mass_, idn_);
  active_num = new_sz;
}

template <typename TypeSystem>
void PointParticles<TypeSystem>::impl_reserve(size_t new_cap) {
  space::reserve_all(new_cap, pos_, vel_, mass_, idn_);
}

template <typename TypeSystem>
size_t PointParticles<TypeSystem>::impl_number() const {
  return active_num;
}

template <typename TypeSystem>
std::ostream &operator<<(std::ostream &os, PointParticles<TypeSystem> const &ps) {
  size_t num = ps.number();
  os << ps.time() << ' ';
  for (size_t i = 0; i < num; ++i) {
    space::display(os, ps.idn()[i], ps.mass()[i], ps.pos().x[i], ps.pos().y[i], ps.pos().z[i], ps.vel().x[i],
                   ps.vel().y[i], ps.vel().z[i]);
  }
  return os;
}
}  // namespace space

#endif  // SPACEHUB_POINT_PARTICLE_H
