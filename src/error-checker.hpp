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
 * @file error-checker.hpp
 *
 * Header file.
 */
#ifndef SPACEHUB_ERROR_CHECKER_HPP
#define SPACEHUB_ERROR_CHECKER_HPP

#include "dev-tools.hpp"

namespace space::ode_iterator {

/*---------------------------------------------------------------------------*\
    Class Error_checker Declaration
\*---------------------------------------------------------------------------*/
/**
 * Abstract class of step error estimator. A class implements(partly/fully) the interfaces of this
 * class via CRTP idiom can be used cross the system as an implementation of the concept `ErrorChecker`. The error estimator
 * provides the interface to estimate the error of an integration system.
 *
 * @tparam Derived The implement class in CRTP idiom.
 */
  template<typename Derived>
  class ErrorChecker {
  public:
    // public methods
    DECLARE_CRTP_READ_ACCESSOR(Derived, auto, atol);

    DECLARE_CRTP_READ_ACCESSOR(Derived, auto, rtol);

    /**
     * @auto_impl
     *
     * The downcast interface of Base class to Derived class.
     * @return Derived
     */
    Derived &derived();

    /**
     * @opt_impl{Absolute error is required.}
     *
     * Set the absolute error of the error estimator.
     *
     * @tparam Scalar Floating point like type type.
     * @param[in] abs_error Absolute error.
     */
    template<typename Scalar>
    void set_atol(Scalar abs_error);

    /**
     * @opt_impl{Relative error is required.}
     *
     * Set the relative error of the error estimator.
     *
     * @tparam Scalar Floating point like type type.
     * @param rel_error Relative error.
     */
    template<typename Scalar>
    void set_rtol(Scalar rel_error);

    /**
     * @must_impl
     *
     * Estimate the error between two results. For example,
     * @f[\mathrm
     *   error = {|y_0 -y_1|/|y_0|}
     * @f]
     *
     * @tparam Array Iterable array like type.
     * @param[in] y0 The first result.
     * @param[in] y1 The second result.
     * @return The estimated error.
     */
    template<typename Array>
    auto error(Array const &y0, Array const &y1) -> typename Array::value_type;

    /**
     * @must_impl
     *
     * Estimate the error between two results with given scale. For example
     * @f[\mathrm
     *   error = {|y_0 -y_1|/scale}
     * @f]
     *
     * @tparam Array Iterable array like type.
     * @param[in] scale The provided scale for the results.
     * @param[in] y0 The first result.
     * @param[in] y1 The second result.
     * @return
     */
    template<typename Array>
    auto error(Array const &scale, Array const &y0, Array const &y1) -> typename Array::value_type;

  private:
    /**
     * Construct a new ErrorChecker object
     */
    ErrorChecker() = default;

    friend Derived;
  };

/*---------------------------------------------------------------------------*\
    Class Particles Implementation
\*---------------------------------------------------------------------------*/
  template<typename Derived>
  Derived &ErrorChecker<Derived>::derived() {
    return static_cast<Derived &>(*this);
  }

  template<typename Derived>
  template<typename T>
  void ErrorChecker<Derived>::set_atol(T atol) {
    static_cast<Derived *>(this)->impl_set_atol(atol);
  }

  template<typename Derived>
  template<typename T>
  void ErrorChecker<Derived>::set_rtol(T rtol) {
    static_cast<Derived *>(this)->impl_set_rtol(rtol);
  }

  template<typename Derived>
  template<typename Array>
  auto ErrorChecker<Derived>::error(Array const &y0, Array const &y1) -> typename Array::value_type {
    return static_cast<Derived *>(this)->impl_error(y0, y1);
  }

  template<typename Derived>
  template<typename Array>
  auto
  ErrorChecker<Derived>::error(Array const &scale, Array const &y0, Array const &y1) -> typename Array::value_type {
    return static_cast<Derived *>(this)->impl_error(scale, y0, y1);
  }

}
#endif //SPACEHUB_ERROR_CHECKER_HPP
