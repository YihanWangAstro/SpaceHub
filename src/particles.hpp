
#ifndef SPACEHUB_PARTICLES_HPP
#define SPACEHUB_PARTICLES_HPP

#include "dev-tools.hpp"

namespace space {

/*---------------------------------------------------------------------------*\
    Class Particles Declaration
\*---------------------------------------------------------------------------*/
/**
 * @brief CRTP base class of a 'Structure of Array' kind particle set.
 *
 * @tparam Derived
 */
  template<typename Derived>
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
     * @brief The particle number of this set.
     *
     * @return size_t
     */
    [[nodiscard]] size_t number() const;

    /**
     *
     * @return
     */
    [[nodiscard]] size_t capacity() const;

    /**
     * @brief Reserve(allocate) space for creating particles
     *
     * @param new_cap
     */
    void reserve(size_t new_cap);

    /**
     * @brief Change the particle number of the set.
     *
     * @param new_sz
     */
    void resize(size_t new_sz);

    /**
     *
     */
    void clear();

    /**
     *
     * @param particle
     */
    template<typename Particle>
    void emplace_back(Particle const &new_particle);

  private:
    /**
     * @brief Construct a new Particles object
     *
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
