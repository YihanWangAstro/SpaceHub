//
// Created by yihan on 2/24/19.
//

#ifndef SPACEHUB_PARTICLE_SYSTEM_HPP
#define SPACEHUB_PARTICLE_SYSTEM_HPP

#include "dev-tools.hpp"

namespace space {
/*---------------------------------------------------------------------------*\
    Class ParticleSystem Declaration
\*---------------------------------------------------------------------------*/
/**
 *
 * @tparam Derived
 */
  template<typename Derived>
  class ParticleSystem {
  public:
    // public methods
    /**
     *
     * @return
     */
    DECLARE_CRTP_ACCESSOR(Derived, auto, mass);

    /**
     *
     * @return
     */
    DECLARE_CRTP_ACCESSOR(Derived, auto, idn);

    /**
     *
     * @return
     */
    DECLARE_CRTP_ACCESSOR(Derived, auto, time);

    /**
     *
     * @return
     */
    DECLARE_CRTP_ACCESSOR(Derived, auto, pos);

    /**
     *
     * @return
     */
    DECLARE_CRTP_ACCESSOR(Derived, auto, vel);

    /**
     *
     * @return
     */
    [[nodiscard]] size_t number() const;

    /**
     *
     * @tparam Scalar
     * @param dt
     */
    template<typename Scalar>
    void advance_time(Scalar dt);

    /**
     *
     * @tparam Coord
     * @tparam Scalar
     * @param step_size
     * @param velocity
     */
    template<typename Coord, typename Scalar>
    void advance_pos(Scalar step_size, Coord const &velocity);

    /**
     *
     * @tparam Coord
     * @tparam Scalar
     * @param step_size
     * @param acceleration
     */
    template<typename Coord, typename Scalar>
    void advance_vel(Scalar step_size, Coord const &acceleration);

    /**
     *
     * @tparam Coord
     * @param acceleration
     */
    template<typename Coord>
    void evaluate_acc(Coord &acceleration) const;

    /**
     *
     * @tparam Scalar
     * @param step_size
     */
    template<typename Scalar>
    void drift(Scalar step_size);

    /**
     *
     * @tparam Scalar
     * @param step_size
     */
    template<typename Scalar>
    void kick(Scalar step_size);

    /**
     *
     */
    void pre_iter_process();

    /**
     *
     */
    void post_iter_process();

    /**
     *
     * @tparam STL
     * @param stl_ranges
     */
    template<typename STL>
    void to_linear_container(STL &stl_ranges);

    /**
     *
     * @tparam STL
     * @param stl_ranges
     */
    template<typename STL>
    void load_from_linear_container(STL const &stl_ranges);

  private:
    // constructors
    ParticleSystem() = default;

    friend Derived;

    // private methods
    void impl_pre_iter_process() {}  // default implementation

    void impl_post_iter_process() {}  // default implementation
  };

/*---------------------------------------------------------------------------*\
    Class ParticleSystem Implementation
\*---------------------------------------------------------------------------*/
  template<typename Derived>
  size_t ParticleSystem<Derived>::number() const {
    return static_cast<Derived const *>(this)->impl_number();
  }

  template<typename Derived>
  template<typename Scalar>
  void ParticleSystem<Derived>::advance_time(Scalar dt) {
    static_cast<Derived *>(this)->impl_advance_time(dt);
  }

  template<typename Derived>
  template<typename Coord, typename Scalar>
  void ParticleSystem<Derived>::advance_pos(Scalar step_size, const Coord &velocity) {
    static_cast<Derived *>(this)->impl_advance_pos(step_size, velocity);
  }

  template<typename Derived>
  template<typename Coord, typename Scalar>
  void ParticleSystem<Derived>::advance_vel(Scalar step_size, const Coord &acceleration) {
    static_cast<Derived *>(this)->impl_advance_vel(step_size, acceleration);
  }

  template<typename Derived>
  template<typename Coord>
  void ParticleSystem<Derived>::evaluate_acc(Coord &acceleration) const {
    static_cast<Derived const *>(this)->impl_evaluate_acc(acceleration);
  }

  template<typename Derived>
  template<typename Scalar>
  void ParticleSystem<Derived>::drift(Scalar step_size) {
    static_cast<Derived *>(this)->impl_drift(step_size);
  }

  template<typename Derived>
  template<typename Scalar>
  void ParticleSystem<Derived>::kick(Scalar step_size) {
    static_cast<Derived *>(this)->impl_kick(step_size);
  }

  template<typename Derived>
  void ParticleSystem<Derived>::pre_iter_process() {
    static_cast<Derived *>(this)->impl_pre_iter_process();
  }

  template<typename Derived>
  void ParticleSystem<Derived>::post_iter_process() {
    static_cast<Derived *>(this)->impl_post_iter_process();
  }

  template<typename Derived>
  template<typename STL>
  void ParticleSystem<Derived>::to_linear_container(STL &stl_ranges) {
    static_assert(is_container_v<STL>, "Only STL-like container can be used");
    static_cast<Derived *>(this)->impl_to_linear_container(stl_ranges);
  }

  template<typename Derived>
  template<typename STL>
  void ParticleSystem<Derived>::load_from_linear_container(const STL &stl_ranges) {
    static_assert(is_container_v<STL>, "Only STL-like container can be used");
    static_cast<Derived *>(this)->impl_load_from_linear_container(stl_ranges);
  }

  template<typename Derived>
  std::ostream &operator<<(std::ostream &os, ParticleSystem<Derived> const &ps) {
    os << static_cast<Derived const &>(ps);
    return os;
  }

  template<typename Derived>
  std::istream &operator>>(std::istream &is, ParticleSystem<Derived> &ps) {
    is >> static_cast<Derived &>(ps);
    return is;
  }

/*---------------------------------------------------------------------------*\
    Help functions and tools
\*---------------------------------------------------------------------------*/
  template<typename T>
  constexpr bool is_particle_system_v = std::is_base_of_v<ParticleSystem<T>, T>;
}  // namespace space

#endif  // SPACEHUB_PARTICLE_SYSTEM_HPP
