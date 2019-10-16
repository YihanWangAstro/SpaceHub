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
 * @file worst-offender.hpp
 *
 * Header file.
 */
#ifndef SPACEHUB_WORST_OFFENDER_HPP
#define SPACEHUB_WORST_OFFENDER_HPP

#include "error-checker.hpp"


namespace space::ode_iterator {

  /**
   *
   * @tparam T
   */
  template<typename T>
  class WorstOffender : public ErrorChecker<WorstOffender<T>> {
  public:
    //Type member
    using Base = ErrorChecker<WorstOffender<T>>;

    using Scalar = T;

    using value_type = T;

    // Constructors
    SPACEHUB_MAKE_CONSTRUCTORS(WorstOffender, default, default, default, default, default);

    WorstOffender(Scalar atol, Scalar rtol) : atol_{atol}, rtol_{rtol} {}

    CRTP_IMPL :
    // CRTP implementation

    SPACEHUB_READ_ACCESSOR(auto, impl_atol, atol_);

    SPACEHUB_READ_ACCESSOR(auto, impl_rtol, rtol_);

    void impl_set_atol(Scalar);

    void impl_set_rtol(Scalar);

    template<typename Array>
    auto impl_error(Array const &y0, Array const &y1) -> typename Array::value_type;

    template<typename Array>
    auto impl_error(Array const &scale, Array const &y0, Array const &y1) -> typename Array::value_type;

  private:
    Scalar atol_{1e-12};

    Scalar rtol_{1e-12};
  };


  template<typename T>
  void WorstOffender<T>::impl_set_atol(Scalar error) {
    atol_ = error;
  }


  template<typename T>
  void WorstOffender<T>::impl_set_rtol(Scalar error) {
    rtol_ = error;
  }

  template<typename T>
  template<typename Array>
  auto WorstOffender<T>::impl_error(const Array &y0, const Array &y1) -> typename Array::value_type {
    size_t const size = y0.size();
    Scalar max_err = 0;
    for (size_t i = 0; i < size; ++i) {
      max_err = std::max(max_err, fabs(y0[i] - y1[i]) / (atol_ + std::max(fabs(y0[i]), fabs(y1[i])) * rtol_));
    }
    return max_err;
  }

  template<typename T>
  template<typename Array>
  auto
  WorstOffender<T>::impl_error(const Array &scale, const Array &y0, const Array &y1) -> typename Array::value_type {
    size_t const size = scale.size();
    Scalar max_err = 0;
    for (size_t i = 0; i < size; ++i) {
      max_err = std::max(max_err, fabs(y0[i] - y1[i]) / (atol_ + fabs(scale[i]) * rtol_));
    }
    return max_err;
  }
}
#endif //SPACEHUB_WORST_OFFENDER_HPP
