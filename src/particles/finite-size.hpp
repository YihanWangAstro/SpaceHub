/*---------------------------------------------------------------------------*\
        .-''''-.         |
       /        \        |
      /_        _\       |  SpaceHub: The Open Source N-body Toolkit
     // \  <>  / \\      |
     |\__\    /__/|      |  Website:  https://yihanwangastro.github.io/SpaceHub/
      \    ||    /       |
        \  __  /         |  Copyright (C) 2019 Yihan Wang
         '.__.'          |
---------------------------------------------------------------------
License
    This file is part of SpaceHub.
    SpaceHub is free software: you can redistribute it and/or modify it under
    the terms of the MIT License. SpaceHub is distributed in the hope that it
    will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the MIT License
    for more details. You should have received a copy of the MIT License along
    with SpaceHub.
\*---------------------------------------------------------------------------*/
/**
 * @file finite-size.hpp
 *
 * Header file.
 */
#ifndef SPACEHUB_FINITE_SIZE_HPP
#define SPACEHUB_FINITE_SIZE_HPP

#include "../IO.hpp"
#include "point-particles.hpp"

namespace space::particle_set {

/*---------------------------------------------------------------------------*\
    Class SizeParticles Declaration
\*---------------------------------------------------------------------------*/
/**
 * @brief Finite size particle.
 *
 * @tparam TypeSystem
 */
template <typename TypeSystem>
class SizeParticles : public Particles<SizeParticles<TypeSystem>> {
 public:
  // Type members
  SPACEHUB_USING_TYPE_SYSTEM_OF(TypeSystem);

  using Base = Particles<SizeParticles<TypeSystem>>;

  /*---------------------------------------------------------------------------*\
  Sub-Class Particle Declaration and Implementation
  \*---------------------------------------------------------------------------*/
  /**
   * @brief Embedded Particle of SOA particles `SizeParticles`.
   *
   */
  struct Particle : public PointParticle<Scalar> {
   public:
    // Type members
    using Scalar = typename TypeSystem::Scalar;

    using Vector = Vec3<Scalar>;

    // Constructors
    SPACEHUB_MAKE_CONSTRUCTORS(Particle, default, default, default, default, default);

    /**
     * @brief Construct a new Particle object
     *
     * @param m
     * @param r
     * @param p
     * @param v
     */
    Particle(Scalar m, Scalar r, Vector p, Vector v) : PointParticle<Scalar>{m, p, v}, radius{r} {}

    /**
     * @brief Construct a new Particle object
     *
     * @param m
     * @param r
     * @param px
     * @param py
     * @param pz
     * @param vx
     * @param vy
     * @param vz
     */
    Particle(Scalar m, Scalar r, Scalar px = 0, Scalar py = 0, Scalar pz = 0, Scalar vx = 0, Scalar vy = 0,
             Scalar vz = 0)
        : PointParticle<Scalar>{m, px, py, pz, vx, vy, vz}, radius{r} {}

    friend std::ostream &operator<<(std::ostream &os, Particle const &particle) {
      space::print_csv(os, particle.mass, particle.radius, particle.pos, particle.vel);
      return os;
    }

    friend std::istream &operator>>(std::istream &is, Particle &particle) {
      space::input(is, particle.mass, particle.radius, particle.pos, particle.vel);
      return is;
    }

    // Public members
    Scalar radius;
  };

  // Constructors
  SPACEHUB_MAKE_CONSTRUCTORS(SizeParticles, default, default, default, default, default);

  template <typename STL>
  SizeParticles(Scalar time, STL const &particles_set);

  // Public methods
  SPACEHUB_STD_ACCESSOR(auto, radius, radius_);

  CRTP_IMPL :
      // CRTP implementation

      SPACEHUB_STD_ACCESSOR(auto, impl_mass, mass_);

  SPACEHUB_STD_ACCESSOR(auto, impl_idn, idn_);

  SPACEHUB_STD_ACCESSOR(auto, impl_time, time_);

  SPACEHUB_STD_ACCESSOR(auto, impl_pos, pos_);

  SPACEHUB_STD_ACCESSOR(auto, impl_vel, vel_);

  void impl_resize(size_t new_sz);

  void impl_reserve(size_t new_cap);

  void impl_emplace_back(Particle const &new_particle);

  [[nodiscard]] size_t impl_number() const;

  [[nodiscard]] size_t impl_capacity() const;

  void impl_clear();

 private:
  // Private members
  Coord pos_;

  Coord vel_;

  ScalarArray mass_;

  ScalarArray radius_;

  IdxArray idn_;

  Scalar time_;

  size_t active_num_{0};
};

/*---------------------------------------------------------------------------*\
    Class SizeParticles Implementation
\*---------------------------------------------------------------------------*/
template <typename TypeSystem>
template <typename STL>
SizeParticles<TypeSystem>::SizeParticles(Scalar time, const STL &particles_set) {
  static_assert(is_ranges_v<STL>, "Only STL-like container can be used");
  SPACEHUB_PARTICLE_TYPE_CHECK(STL, Particle);

  size_t input_num = particles_set.size();
  this->reserve(input_num);
  size_t id = 0;
  for (auto &p : particles_set) {
    pos_.emplace_back(p.pos);
    vel_.emplace_back(p.vel);
    mass_.emplace_back(p.mass);
    radius_.emplace_back(p.radius);
    idn_.emplace_back(id++);
  }
  time_ = time;
  active_num_ = input_num;
}

template <typename TypeSystem>
size_t SizeParticles<TypeSystem>::impl_number() const {
  return active_num_;
}

template <typename TypeSystem>
size_t SizeParticles<TypeSystem>::impl_capacity() const {
  return idn_.capacity();
}

template <typename TypeSystem>
void SizeParticles<TypeSystem>::impl_reserve(size_t new_cap) {
  space::reserve_all(new_cap, pos_, vel_, mass_, radius_, idn_);
}

template <typename TypeSystem>
void SizeParticles<TypeSystem>::impl_clear() {
  space::clear_all(pos_, vel_, mass_, radius_, idn_);
  active_num_ = 0;
}

template <typename TypeSystem>
void SizeParticles<TypeSystem>::impl_resize(size_t new_sz) {
  space::resize_all(new_sz, pos_, vel_, mass_, radius_, idn_);
  active_num_ = new_sz;
}

template <typename TypeSystem>
void SizeParticles<TypeSystem>::impl_emplace_back(typename SizeParticles<TypeSystem>::Particle const &new_particle) {
  pos_.emplace_back(new_particle.pos);
  vel_.emplace_back(new_particle.vel);
  mass_.emplace_back(new_particle.mass);
  radius_.emplace_back(new_particle.radius);
  idn_.emplace_back(this->number());
  active_num_++;
}

template <typename TypeSystem>
std::ostream &operator<<(std::ostream &os, SizeParticles<TypeSystem> const &ps) {
  size_t num = ps.number();
  os << ps.time();
  for (size_t i = 0; i < num; ++i) {
    space::print_csv(os, "", ps.idn()[i], ps.mass()[i], ps.radius()[i], ps.pos().x[i], ps.pos().y[i], ps.pos().z[i],
                     ps.vel().x[i], ps.vel().y[i], ps.vel().z[i]);
  }
  return os;
}
}  // namespace space::particle_set
#endif  // SPACEHUB_FINITE_SIZE_HPP
