
#ifndef SPACEHUB_PARTICLES_HPP
#define SPACEHUB_PARTICLES_HPP

#include "dev-tools.hpp"

namespace space {

/*---------------------------------------------------------------------------*\
    Class Particles Declaration
\*---------------------------------------------------------------------------*/
/**
 * Abstract class of a 'Structure of Array(SoA)' particle set. A class implements(partly/fully) the interfaces of this
 * class via CRTP idiom can be used cross the system as an implementation of the concept `Particles`.
 *
 * @tparam Derived The implement class in CRTP idiom.
 */
  template<typename Derived>
  class Particles {
  public:
    // public methods
    DECLARE_CRTP_ACCESSOR(Derived, auto, idn);

    DECLARE_CRTP_ACCESSOR(Derived, auto, mass);

    DECLARE_CRTP_ACCESSOR(Derived, auto, pos);

    DECLARE_CRTP_ACCESSOR(Derived, auto, time);

    DECLARE_CRTP_ACCESSOR(Derived, auto, vel);

    /**
     * @glabel{auto impl}
     *
     * The downcast interface of Base class to Derived class.
     * @return Derived
     */
    Derived &derived();

    /**
     * @must_impl
     *
     * Get the (active) particle number of the SoA particle set.
     * @return size_t the number of the particles.
     */
    [[nodiscard]] size_t number() const;

    /**
     * @opt_impl{The particle set is resizable.}
     *
     * Get the capacity of the `Container` for the each component(i.e idn, pos, vel) of the SoA particle set.
     * @return size_t the capacity of the the SoA particle set.
     * @note implementation should keep the consistence of the capacity of all components.
     */
    [[nodiscard]] size_t capacity() const;

    /**
     * @opt_impl{The particle set is resizable.}
     *
     * Reserve(allocate) space for each component(i.e idn, pos, vel) of the SoA particle set.
     *
     * @param[in] new_cap The new capacity size(in particle number) of the SoA particle set.
     */
    void reserve(size_t new_cap);

    /**
     * @opt_impl{The particle set is resizable.}
     *
     * Change the particle number of the set.
     *
     * @param[in] new_sz The new size(in particle number) of the SoA particle set.
     */
    void resize(size_t new_sz);

    /**
     * @must_impl
     *
     * Clear all the particles in the set.
     * @note After the clear, the number of the particle set should be 0. i.e number() should return 0;
     */
    void clear();

    /**
     * @must_impl
     *
     * Add a single particle to the SoA particle set.
     * @tparam Particle The generic particle Type of the added particle.
     * @param[in] new_particle The single particle that is going to be added.
     */
    template<typename Particle>
    void emplace_back(Particle const &new_particle);

  private:
    /**
     * Construct a new Particles object. Set to be private to avoid outside access. Only implement class can access.
     */
    Particles() = default;

    friend Derived;
  };

/*---------------------------------------------------------------------------*\
    Class Particles Implementation
\*---------------------------------------------------------------------------*/
  template<typename Derived>
  Derived &Particles<Derived>::derived() {
    return static_cast<Derived &>(*this);
  }

  template<typename Derived>
  size_t Particles<Derived>::number() const {
    return static_cast<Derived const *>(this)->impl_number();
  }

  template<typename Derived>
  size_t Particles<Derived>::capacity() const {
    return static_cast<Derived const *>(this)->impl_capacity();
  }

  template<typename Derived>
  void Particles<Derived>::reserve(size_t new_cap) {
    static_cast<Derived *>(this)->impl_reserve(new_cap);
  }

  template<typename Derived>
  void Particles<Derived>::resize(size_t new_sz) {
    static_cast<Derived *>(this)->impl_resize(new_sz);
  }

  template<typename Derived>
  void Particles<Derived>::clear() {
    static_cast<Derived *>(this)->impl_clear();
  }

  template<typename Derived>
  template<typename Particle>
  void Particles<Derived>::emplace_back(Particle const &particle) {
    static_cast<Derived *>(this)->impl_emplace_back(particle);
  }

/*---------------------------------------------------------------------------*\
    Help functions and tools
\*---------------------------------------------------------------------------*/

  template<typename T>
  constexpr bool is_soa_particles_v = std::is_base_of_v<Particles<T>, T>;

#define SPACEHUB_PARTICLE_TYPE_CHECK(CTR, VAL)                    \
  static_assert(std::is_base_of_v<typename CTR::value_type, VAL>, \
                "Class can only be initialized by containers with its internal 'Particle' type!");

}  // namespace space
#endif
