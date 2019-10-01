//
// Created by yihan on 1/3/19.
//

#ifndef SPACEHUB_EOM_H
#define SPACEHUB_EOM_H

#include "core-computation.hpp"
#include "type-class.hpp"

namespace space {

/*---------------------------------------------------------------------------*\
    Class Interactions Declaration
\*---------------------------------------------------------------------------*/
template <typename Derived>
class Interactions {
 private:
  // special macros to generate compile time method check
  CREATE_METHOD_CHECK(impl_eval_extra_vel_indep_acc);

  CREATE_METHOD_CHECK(impl_eval_extra_vel_dep_acc);

 public:
  // public static members
  /**
   *
   */
  static constexpr bool ext_vel_dep{HAS_METHOD(Derived, impl_eval_extra_vel_dep_acc)};

  /**
   *
   */
  static constexpr bool ext_vel_indep{HAS_METHOD(Derived, impl_eval_extra_vel_indep_acc)};

  // public method
  /**
   *
   * @tparam Particles
   * @param partc
   * @param acc
   */
  template <typename Particles>
  void eval_acc(Particles const &partc, typename Particles::Coord &acc);

  /**
   *
   * @tparam Particles
   * @param partc
   * @param acc
   */
  template <typename Particles>
  void eval_extra_vel_dep_acc(Particles const &partc, typename Particles::Coord &acc);

  /**
   *
   * @tparam Particles
   * @param partc
   * @param acc
   */
  template <typename Particles>
  void eval_extra_vel_indep_acc(Particles const &partc, typename Particles::Coord &acc);

  /**
   *
   * @tparam Particles
   * @param partc
   * @param acc
   */
  template <typename Particles>
  void eval_newtonian_acc(Particles const &partc, typename Particles::Coord &acc);

 private:
  // constructors
  Interactions() = default;

  friend Derived;
};

/*---------------------------------------------------------------------------*\
    Class Interactions Implementation
\*---------------------------------------------------------------------------*/
template <typename Derived>
template <typename Particles>
void Interactions<Derived>::eval_acc(const Particles &partc, typename Particles::Coord &acc) {
  static_cast<Derived *>(this)->impl_eval_acc(partc, acc);
}

template <typename Derived>
template <typename Particles>
void Interactions<Derived>::eval_extra_vel_dep_acc(const Particles &partc, typename Particles::Coord &acc) {
  static_cast<Derived *>(this)->impl_eval_extra_vel_dep_acc(partc, acc);
}

template <typename Derived>
template <typename Particles>
void Interactions<Derived>::eval_extra_vel_indep_acc(const Particles &partc, typename Particles::Coord &acc) {
  static_cast<Derived *>(this)->impl_eval_extra_vel_indep_acc(partc, acc);
}

template <typename Derived>
template <typename Particles>
void Interactions<Derived>::eval_newtonian_acc(const Particles &partc, typename Particles::Coord &acc) {
  static_cast<Derived *>(this)->impl_eval_newtonian_acc(partc, acc);
}
/*---------------------------------------------------------------------------*\
    Help functions and tools
\*---------------------------------------------------------------------------*/
template <typename T>
constexpr bool is_interactions_v = std::is_base_of_v<Interactions<T>, T>;
}  // namespace space

#endif  // SPACEHUB_EOM_H