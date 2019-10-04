
#ifndef PARTICLES_H
#define PARTICLES_H

#include "dev-tools.hpp"
#include "type-class.hpp"
#include "vector/vector3.h"

namespace space {

/*---------------------------------------------------------------------------*\
    Class PointParticle Declaration
\*---------------------------------------------------------------------------*/
template <typename Real>
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
   * @param mass
   * @param position
   * @param velocity
   */
  explicit PointParticle(Scalar mass, Vector position, Vector velocity);

  /**
   * @brief Construct a new Point Particle object
   *
   * @param mass
   * @param px
   * @param py
   * @param pz
   * @param vx
   * @param vy
   * @param vz
   */
  explicit PointParticle(Scalar mass, Scalar px = 0, Scalar py = 0, Scalar pz = 0, Scalar vx = 0, Scalar vy = 0,
                         Scalar vz = 0);

  // public members
  Vector pos;

  Vector vel;

  Scalar mass;
};

/*---------------------------------------------------------------------------*\
    Class Particles Declaration
\*---------------------------------------------------------------------------*/
template <typename Derived>
class Particles {
 public:
  // public methods
  /**
   *
   * @return
   */
  DECLARE_CRTP_ACCESSOR(Derived, auto, idn);

  /**
   *
   * @return
   */
  DECLARE_CRTP_ACCESSOR(Derived, auto, mass);

  /**
   *
   * @return
   */
  DECLARE_CRTP_ACCESSOR(Derived, auto, pos);

  /**
   *
   * @return
   */
  DECLARE_CRTP_ACCESSOR(Derived, auto, time);

  /**
   *
   * @return
   */
  DECLARE_CRTP_ACCESSOR(Derived, auto, vel);

  /**
   *
   * @return
   */
  Derived &derived();

  /**
   *
   * @return
   */
  size_t number() const;

  /**
   *
   * @param new_cap
   */
  void reserve(size_t new_cap);

  /**
   *
   * @param new_sz
   */
  void resize(size_t new_sz);

 private:
  /**
   * @brief Construct a new Particles object
   *
   */
  Particles() = default;

  friend Derived;
};

/*---------------------------------------------------------------------------*\
    Class PointParticle Implementation
\*---------------------------------------------------------------------------*/
template <typename Real>
PointParticle<Real>::PointParticle(Scalar m, PointParticle::Vector p, PointParticle::Vector v)
    : pos(p), vel(v), mass(m) {}

template <typename Real>
PointParticle<Real>::PointParticle(Scalar m, Scalar px, Scalar py, Scalar pz, Scalar vx, Scalar vy, Scalar vz)
    : pos(px, py, pz), vel(vx, vy, vz), mass(m) {}

template <typename Real>
std::ostream &operator<<(std::ostream &os, PointParticle<Real> const &particle) {
  space::display(os, particle.mass, particle.pos, particle.vel);
  return os;
}

template <typename Real>
std::istream &operator>>(std::istream &is, PointParticle<Real> &particle) {
  space::input(is, particle.mass, particle.pos, particle.vel);
  return is;
}

/*---------------------------------------------------------------------------*\
    Class Particles Implementation
\*---------------------------------------------------------------------------*/
template <typename Derived>
Derived &Particles<Derived>::derived() {
  return static_cast<Derived &>(*this);
}

template <typename Derived>
size_t Particles<Derived>::number() const {
  return static_cast<Derived const *>(this)->impl_number();
}

template <typename Derived>
void Particles<Derived>::reserve(size_t new_cap) {
  static_cast<Derived *>(this)->impl_reserve(new_cap);
}

template <typename Derived>
void Particles<Derived>::resize(size_t new_sz) {
  static_cast<Derived *>(this)->impl_resize(new_sz);
}

/*---------------------------------------------------------------------------*\
    Help functions and tools
\*---------------------------------------------------------------------------*/
template <typename T>
constexpr bool is_soa_particles_v = std::is_base_of_v<Particles<T>, T>;

#define SPACEHUB_PARTICLE_TYPE_CHECK(CTR, VAL)                    \
  static_assert(std::is_base_of_v<typename CTR::value_type, VAL>, \
                "Class can only be initialized by containers with its internal 'Particle' type!");

}  // namespace space
#endif
