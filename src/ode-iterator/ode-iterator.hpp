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
 * @file ode-iterator.hpp
 *
 * Header file.
 */
#pragma once

/**
 * @namespace space::ode_iterator
 * space name for ode-iterator
 */
namespace space::ode_iterator {

/*---------------------------------------------------------------------------*\
    Class OdeIterator Declaration
\*---------------------------------------------------------------------------*/
/**
 * Abstract class of ode-iterator. A class implements(partly/fully) the interfaces of this
 * class via CRTP idiom can be used cross the system as an implementation of the concept `OdeIterator`.
 * The OdeIterator is used to iterate the ParticleSystem by steps.
 *
 * @tparam Derived The implement class in CRTP idiom.
 */
template <typename Derived>
class OdeIterator {
 public:
  // Public methods
  /**
   * @must_impl
   *
   * Evolve a particle system with given step size as a trial. The method should return the proper
   * step size for next iteration.
   *
   * @tparam T Any implementation of concept `ParticleSystem`.
   * @param[in,out] particles The particles that need evolution.
   * @param[in] macro_step_size Step size.
   * @return The step size for next iteration.
   */
  template <typename T>
  auto iterate(T &particles, typename T::Scalar macro_step_size) -> typename T::Scalar;

  /**
   * @opt_impl{Absolute error tolerance is required.}
   *
   * Set the absolute error tolerance of the iterator.
   *
   * @tparam Scalar Floating point like scalar.
   * @param[in] atol Absolute error tolerance.
   */
  template <typename Scalar>
  void set_atol(Scalar atol);

  /**
   * @opt_impl{Relative error tolerance is required.}
   *
   * Set the relative error tolerance of the iterator.
   *
   * @tparam Scalar Floating point like scalar.
   * @param[in] rtol Relative error tolerance.
   */
  template <typename Scalar>
  void set_rtol(Scalar rtol);

 private:
  // Constructor
  OdeIterator() = default;

  friend Derived;
};

/*---------------------------------------------------------------------------*\
    Class OdeIterator Implementation
\*---------------------------------------------------------------------------*/
template <typename Derived>
template <typename T>
auto OdeIterator<Derived>::iterate(T &particles, typename T::Scalar macro_step_size) -> typename T::Scalar {
  return static_cast<Derived *>(this)->impl_iterate(particles, macro_step_size);
}

template <typename Derived>
template <typename T>
void OdeIterator<Derived>::set_atol(T atol) {
  return static_cast<Derived *>(this)->impl_set_atol(atol);
}

template <typename Derived>
template <typename T>
void OdeIterator<Derived>::set_rtol(T rtol) {
  return static_cast<Derived *>(this)->impl_set_rtol(rtol);
}

/*---------------------------------------------------------------------------*\
    Help functions and tools
\*---------------------------------------------------------------------------*/
template <typename T>
constexpr bool is_ode_iterator_v = std::is_base_of_v<OdeIterator<T>, T>;
}  // namespace space::ode_iterator
