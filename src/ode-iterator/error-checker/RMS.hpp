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
 * @file RMS.hpp
 *
 * Header file.
 */
#ifndef SPACEHUB_RMS_HPP
#define SPACEHUB_RMS_HPP

#include "error-checker.hpp"

namespace space::ode_iterator {

  /**
   *
   * @tparam T
   */
  template<typename T>
  class RMS : public ErrorChecker<RMS<T>> {
  public:
    //Type member
    using Base = ErrorChecker<RMS<T>>;

    using Scalar = T;

    using value_type = T;

    // Constructors
    SPACEHUB_MAKE_CONSTRUCTORS(RMS, default, default, default, default, default);

    RMS(Scalar atol, Scalar rtol) : atol_{atol}, rtol_{rtol} {}

    CRTP_IMPL :
    // CRTP implementation

    SPACEHUB_READ_ACCESSOR(auto, impl_atol, atol_);

    SPACEHUB_READ_ACCESSOR(auto, impl_rtol, rtol_);

    void impl_set_atol(Scalar);

    void impl_set_rtol(Scalar);

    template<typename Array>
    auto impl_error(Array const &scale, Array const &diff) -> typename Array::value_type;

    template<typename Array>
    auto impl_error(Array const &scale, Array const &y0, Array const &y1) -> typename Array::value_type;

  private:
    Scalar atol_{1e-12};

    Scalar rtol_{1e-12};
  private:
    CREATE_MEMBER_CHECK(x);
    CREATE_MEMBER_CHECK(y);
    CREATE_MEMBER_CHECK(z);

    template<typename Array>
    auto one_dimension_error(Array const &scale, Array const &diff);

    template<typename Array>
    auto one_dimension_error(Array const &scale, Array const &y0, Array const &y1);
  };

  template<typename T>
  void RMS<T>::impl_set_atol(Scalar error) {
    atol_ = error;
  }

  template<typename T>
  void RMS<T>::impl_set_rtol(Scalar error) {
    rtol_ = error;
  }

  template<typename T>
  template<typename Array>
  auto RMS<T>::impl_error(const Array &scale, const Array &diff) -> typename Array::value_type {
    if constexpr (HAS_MEMBER(Array, x) && HAS_MEMBER(Array, y) && HAS_MEMBER(Array, z)){
      return sqrt((one_dimension_error(scale.x, diff.x) + one_dimension_error(scale.y, diff.y) + one_dimension_error(scale.z, diff.z)) / 3);
    } else {
      return sqrt(one_dimension_error(scale, diff));
    }
  }

  template<typename T>
  template<typename Array>
  auto RMS<T>::impl_error(const Array &scale, const Array &y0, const Array &y1) -> typename Array::value_type {
    if constexpr (HAS_MEMBER(Array, x) && HAS_MEMBER(Array, y) && HAS_MEMBER(Array, z)){
      return sqrt((one_dimension_error(scale.x, y0.x, y1.x) + one_dimension_error(scale.y, y0.y, y1.y) + one_dimension_error(scale.z, y0.z, y1.z)) / 3);
    } else {
      return sqrt(one_dimension_error(scale, y0, y1));
    }
  }

  template<typename T>
  template<typename Array>
  auto RMS<T>::one_dimension_error(const Array &scale, const Array &diff) {
    size_t const size = scale.size();
    Scalar error = 0;
    for (size_t i = 0; i < size; ++i) {
      auto r = fabs(diff[i]) / (atol_ + fabs(scale[i]) * rtol_);
      error += r * r;
    }
    return error / size;
  }

  template<typename T>
  template<typename Array>
  auto RMS<T>::one_dimension_error(const Array &scale, const Array &y0, const Array &y1) {
    size_t const size = scale.size();
    Scalar error = 0;
    for (size_t i = 0; i < size; ++i) {
      auto r = fabs(y0[i] - y1[i]) / (atol_ + fabs(scale[i]) * rtol_);
      error += r * r;
    }
    return error / size;
  }
}
#endif //SPACEHUB_RMS_HPP
