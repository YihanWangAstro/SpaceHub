//
// Created by yihan on 4/3/19.
//

#ifndef SPACEHUB_POINT_PARTICLE_H
#define SPACEHUB_POINT_PARTICLE_H

#include "../particles.hpp"

namespace space {

/*---------------------------------------------------------------------------*\
Class PointParticle Declaration
\*---------------------------------------------------------------------------*/
  template<typename Real>
  struct PointParticle {
  public:
    // type members
    using Scalar = Real;

    using Vector = Vec3<Scalar>;

    /**
     * @brief Construct a new Point Particle object
     *
     */
    PointParticle() = default;

    /**
     * @brief Construct a new Point Particle object
     *
     * @param mass The mass of the particle
     * @param position The 3d vector position of the particle
     * @param velocity The 3d vector velocity of the particle
     */
    explicit PointParticle(Scalar mass, Vector position, Vector velocity);

    /**
     * @brief Construct a new Point Particle object
     *
     * @param mass The mass fof the particle
     * @param px The x-component of the position vector
     * @param py The y-component of the position vector
     * @param pz The z-component of the position vector
     * @param vx The x-component of the velocity vector
     * @param vy The y-component of the velocity vector
     * @param vz The z-component of the velocity vector
     */
    explicit PointParticle(Scalar mass, Scalar px = 0, Scalar py = 0, Scalar pz = 0, Scalar vx = 0, Scalar vy = 0,
                           Scalar vz = 0);

    // public members
    Vector pos;

    Vector vel;

    Scalar mass;
  };

/*---------------------------------------------------------------------------*\
    Class PointParticles Declaration
\*---------------------------------------------------------------------------*/
  template<typename TypeSystem>
  class PointParticles : public Particles<PointParticles<TypeSystem>> {
  public:
    // Type members
    SPACEHUB_USING_TYPE_SYSTEM_OF(TypeSystem);

    using Base = Particles<PointParticles<TypeSystem>>;

    using Particle = PointParticle<Scalar>;

    // Constructors
    SPACEHUB_MAKE_CONSTRUCTORS(PointParticles, default, default, default, default, default);

    template<typename STL>
    PointParticles(Scalar t, STL const &particle_set);

    CRTP_impl :
    // CRTP implementation

    SPACEHUB_STD_ACCESSOR(auto, impl_mass, mass_);

    SPACEHUB_STD_ACCESSOR(auto, impl_idn, idn_);

    SPACEHUB_STD_ACCESSOR(auto, impl_time, time_);

    SPACEHUB_STD_ACCESSOR(auto, impl_pos, pos_);

    SPACEHUB_STD_ACCESSOR(auto, impl_vel, vel_);

    void impl_resize(size_t new_sz);

    void impl_reserve(size_t new_cap);

    void impl_emplace_back(Particle const &new_particle);

    size_t impl_number() const;

    size_t impl_capacity() const;

    void impl_clear();

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
    Class PointParticle Implementation
\*---------------------------------------------------------------------------*/
  template<typename Real>
  PointParticle<Real>::PointParticle(Scalar m, PointParticle::Vector p, PointParticle::Vector v)
          : pos(p), vel(v), mass(m) {}

  template<typename Real>
  PointParticle<Real>::PointParticle(Scalar m, Scalar px, Scalar py, Scalar pz, Scalar vx, Scalar vy, Scalar vz)
          : pos(px, py, pz), vel(vx, vy, vz), mass(m) {}

  template<typename Real>
  std::ostream &operator<<(std::ostream &os, PointParticle<Real> const &particle) {
    space::print_csv(os, particle.mass, particle.pos, particle.vel);
    return os;
  }

  template<typename Real>
  std::istream &operator>>(std::istream &is, PointParticle<Real> &particle) {
    space::input(is, particle.mass, particle.pos, particle.vel);
    return is;
  }

/*---------------------------------------------------------------------------*\
    Class PointParticles Implementation
\*---------------------------------------------------------------------------*/
  template<typename TypeSystem>
  template<typename STL>
  PointParticles<TypeSystem>::PointParticles(Scalar t, const STL &particle_set) {
    static_assert(is_container_v<STL>, "Only STL-like container can be used");
    SPACEHUB_PARTICLE_TYPE_CHECK(STL, Particle);

    size_t input_num = particle_set.size();
    this->reserve(input_num);
    size_t id = 0;
    for (auto &p : particle_set) {
      pos_.emplace_back(p.pos);
      vel_.emplace_back(p.vel);
      mass_.emplace_back(p.mass);
      idn_.emplace_back(id++);
    }
    time_ = t;
    active_num = input_num;
  }

  template<typename TypeSystem>
  void PointParticles<TypeSystem>::impl_resize(size_t new_sz) {
    space::resize_all(new_sz, pos_, vel_, mass_, idn_);
    active_num = new_sz;
  }

  template<typename TypeSystem>
  void PointParticles<TypeSystem>::impl_reserve(size_t new_cap) {
    space::reserve_all(new_cap, pos_, vel_, mass_, idn_);
  }

  template<typename TypeSystem>
  void PointParticles<TypeSystem>::impl_clear() {
    space::clear_all(pos_, vel_, mass_, idn_);
    active_num = 0;
  }

  template<typename TypeSystem>
  void
  PointParticles<TypeSystem>::impl_emplace_back(typename PointParticles<TypeSystem>::Particle const &new_particle) {
    pos_.emplace_back(new_particle.pos);
    vel_.emplace_back(new_particle.vel);
    mass_.emplace_back(new_particle.mass);
    idn_.emplace_back(this->number());
    active_num++;
  }

  template<typename TypeSystem>
  size_t PointParticles<TypeSystem>::impl_number() const {
    return active_num;
  }

  template<typename TypeSystem>
  size_t PointParticles<TypeSystem>::impl_capacity() const {
    return idn_.capacity();
  }

  template<typename TypeSystem>
  std::ostream &operator<<(std::ostream &os, PointParticles<TypeSystem> const &particle_system) {
    size_t num = particle_system.number();
    os << particle_system.time();
    for (size_t i = 0; i < num; ++i) {
      space::print_csv(os, "", particle_system.idn()[i], particle_system.mass()[i], particle_system.pos().x[i],
                       particle_system.pos().y[i], particle_system.pos().z[i], particle_system.vel().x[i],
                       particle_system.vel().y[i], particle_system.vel().z[i]);
    }
    return os;
  }
}  // namespace space

#endif  // SPACEHUB_POINT_PARTICLE_H
