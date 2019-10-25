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
 * @file interaction.hpp
 *
 * Header file.
 */
#ifndef SPACEHUB_INTERACTION_HPP
#define SPACEHUB_INTERACTION_HPP

#include "../core-computation.hpp"
#include "../type-class.hpp"

namespace space::interactions {

/*---------------------------------------------------------------------------*\
    Class Interactions Declaration
\*---------------------------------------------------------------------------*/
/**
 * Abstract class of interactions. A class implements(partly/fully) the interfaces of this
 * class via CRTP idiom can be used cross the system as interaction of the concept `Interactions`.
 *
 * @tparam Derived The implement class in CRTP idiom.
 */
template <typename Derived>
class Interactions {
 private:
  // special macros to generate compile time method check
  CREATE_METHOD_CHECK(impl_eval_extra_vel_indep_acc);

  CREATE_METHOD_CHECK(impl_eval_extra_vel_dep_acc);

 public:
  // public static members
  /**
   * @auto_impl
   *
   * Variable to check if external velocity dependent forces exist the interaction. If the Derived class implements
   * method impl_eval_extra_vel_dep_acc(). The value is `true`, otherwise this value is `false`.
   */
  static constexpr bool ext_vel_dep{HAS_METHOD(Derived, impl_eval_extra_vel_dep_acc)};

  /**
   * @auto_impl
   *
   * Variable to check if external velocity independent forces exist the interaction. If the Derived class implements
   * method impl_eval_extra_vel_indep_acc(). The value is `true`, otherwise this value is `false`.
   */
  static constexpr bool ext_vel_indep{HAS_METHOD(Derived, impl_eval_extra_vel_indep_acc)};

  // public method
  /**
   * @auto_impl
   *
   * The downcast interface of Base class to Derived class.
   * @return Derived
   */
  Derived &derived();

  /**
   * @must_impl
   *
   * Evaluate the total acceleration of the current state of a given particle system.
   *
   * @tparam Particles Type of the particle system. This type must have method mass(), pos(), vel(), time(); type member
   * `Coord`.
   *
   * @param[in] particles The particle system need to be evaluated.
   * @param[out] acceleration The output of the evaluated acceleration.
   */
  template <typename Particles>
  void eval_acc(Particles const &particles, typename Particles::Coord &acceleration) const;

  /**
   * @opt_impl{external velocity dependent forces exist.}
   *
   * Evaluate the external velocity dependent acceleration of the current state of a given particle system.
   *
   * @tparam Particles Type of the particle system. This type must have method mass(), pos(), vel(); type member
   * `Coord`.
   *
   * @param[in] particles The particle system need to be evaluated.
   * @param[out] acceleration The output of the evaluated acceleration.
   */
  template <typename Particles>
  void eval_extra_vel_dep_acc(Particles const &particles, typename Particles::Coord &acceleration) const;

  /**
   * @opt_impl{external velocity independent forces exist.}
   *
   * Evaluate the external velocity independent acceleration of the current state of a given particle system.
   *
   * @tparam Particles Type of the particle system. This type must have method mass(), pos(); type member `Coord`.
   *
   * @param[in] particles The particle system need to be evaluated.
   * @param[out] acceleration The output of the evaluated acceleration.
   */
  template <typename Particles>
  void eval_extra_vel_indep_acc(Particles const &particles, typename Particles::Coord &acceleration) const;

  /**
   * @must_impl
   *
   * Evaluate the internal newtonian acceleration of the current state of a given particle system.
   *
   * @tparam Particles Type of the particle system. This type must have method mass(), pos(); type member `Coord`.
   *
   * @param[in] particles The particle system need to be evaluated.
   * @param[out] acceleration The output of the evaluated acceleration.
   */
  template <typename Particles>
  void eval_newtonian_acc(Particles const &particles, typename Particles::Coord &acceleration) const;

 private:
  // constructors
  Interactions() = default;

  friend Derived;
};

/*---------------------------------------------------------------------------*\
    Class Interactions Implementation
\*---------------------------------------------------------------------------*/
template <typename Derived>
Derived &Interactions<Derived>::derived() {
  return static_cast<Derived &>(*this);
}

template <typename Derived>
template <typename Particles>
void Interactions<Derived>::eval_acc(const Particles &particles, typename Particles::Coord &acceleration) const {
  static_cast<Derived const *>(this)->impl_eval_acc(particles, acceleration);
}

template <typename Derived>
template <typename Particles>
void Interactions<Derived>::eval_extra_vel_dep_acc(const Particles &particles,
                                                   typename Particles::Coord &acceleration) const {
  static_cast<Derived const *>(this)->impl_eval_extra_vel_dep_acc(particles, acceleration);
}

template <typename Derived>
template <typename Particles>
void Interactions<Derived>::eval_extra_vel_indep_acc(const Particles &particles,
                                                     typename Particles::Coord &acceleration) const {
  static_cast<Derived const *>(this)->impl_eval_extra_vel_indep_acc(particles, acceleration);
}

template <typename Derived>
template <typename Particles>
void Interactions<Derived>::eval_newtonian_acc(const Particles &particles,
                                               typename Particles::Coord &acceleration) const {
  static_cast<Derived const *>(this)->impl_eval_newtonian_acc(particles, acceleration);
}

/*---------------------------------------------------------------------------*\
    Help functions and tools
\*---------------------------------------------------------------------------*/
template <typename T>
constexpr bool is_interactions_v = std::is_base_of_v<Interactions<T>, T>;
}  // namespace space::interactions

#endif  // SPACEHUB_INTERACTION_HPP