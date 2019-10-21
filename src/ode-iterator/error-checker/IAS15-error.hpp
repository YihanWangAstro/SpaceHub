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
 * @file IAS15-error.hpp
 *
 * Header file.
 */
#ifndef SPACEHUB_IAS15_ERROR_HPP
#define SPACEHUB_IAS15_ERROR_HPP
#include "error-checker.hpp"

namespace space::ode_iterator {

  /**
   *
   * @tparam T
   */
  template<typename T>
  class IAS15Error : public ErrorChecker<IAS15Error<T>> {
  public:
    //Type member
    using Base = ErrorChecker<IAS15Error<T>>;

    using Scalar = T;

    using value_type = T;

    // Constructors
    SPACEHUB_MAKE_CONSTRUCTORS(IAS15Error, default, default, default, default, default);

    IAS15Error(Scalar atol, Scalar rtol) : atol_{atol}, rtol_{rtol} {}

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

    CREATE_MEMBER_CHECK(x);
    CREATE_MEMBER_CHECK(y);
    CREATE_MEMBER_CHECK(z);

    template<typename Array>
    auto one_dimension_error(Array const &scale, Array const &diff);

    template<typename Array>
    auto one_dimension_error(Array const &scale, Array const &y0, Array const &y1);

  };

  template<typename T>
  void IAS15Error<T>::impl_set_atol(Scalar error) {
    atol_ = error;
  }

  template<typename T>
  void IAS15Error<T>::impl_set_rtol(Scalar error) {
    rtol_ = error;
  }

  template<typename T>
  template<typename Array>
  auto IAS15Error<T>::impl_error(const Array &scale, const Array &diff) -> typename Array::value_type {
    if constexpr (HAS_MEMBER(Array, x) && HAS_MEMBER(Array, y) && HAS_MEMBER(Array, z)){
      auto [max_diff_x, max_scale_x] = one_dimension_error(scale.x, diff.x);
      auto [max_diff_y, max_scale_y] = one_dimension_error(scale.y, diff.y);
      auto [max_diff_z, max_scale_z] = one_dimension_error(scale.z, diff.z);
      return std::max(std::max(max_diff_x, max_diff_y), max_diff_z) / std::max(std::max(max_scale_x, max_scale_y), max_scale_z);
    } else {
      auto [max_diff, max_scale] = one_dimension_error(scale, diff);
      return max_diff / max_scale;
    }
  }

  template<typename T>
  template<typename Array>
  auto
  IAS15Error<T>::impl_error(const Array &scale, const Array &y0, const Array &y1) -> typename Array::value_type {
    if constexpr (HAS_MEMBER(Array, x) && HAS_MEMBER(Array, y) && HAS_MEMBER(Array, z)){
      auto [max_diff_x, max_scale_x] = one_dimension_error(scale.x, y0.x, y1.x);
      auto [max_diff_y, max_scale_y] = one_dimension_error(scale.y, y0.y, y1.y);
      auto [max_diff_z, max_scale_z] = one_dimension_error(scale.z, y0.z, y1.z);
      return std::max(std::max(max_diff_x, max_diff_y), max_diff_z) / std::max(std::max(max_scale_x, max_scale_y), max_scale_z);
    } else {
      auto [max_diff, max_scale] = one_dimension_error(scale, y0, y1);
      return max_diff / max_scale;
    }
  }

  template<typename T>
  template<typename Array>
  auto IAS15Error<T>::one_dimension_error(const Array &scale, const Array &diff) {
    size_t const size = scale.size();
    Scalar max_diff  = 0;
    Scalar max_scale =0;
    for (size_t i = 0; i < size; ++i) {
      max_diff = std::max(max_diff, static_cast<Scalar>(fabs(diff[i])));
      max_scale= std::max(max_scale, static_cast<Scalar>(atol_ + fabs(scale[i]) * rtol_));
    }
    return std::make_tuple(max_diff, max_scale);
  }

  template<typename T>
  template<typename Array>
  auto
  IAS15Error<T>::one_dimension_error(const Array &scale, const Array &y0, const Array &y1) {
    size_t const size = scale.size();
    Scalar max_diff  = 0;
    Scalar max_scale = 0;
    for (size_t i = 0; i < size; ++i) {
      max_diff = std::max(max_diff, static_cast<Scalar>(fabs(y0[i] - y1[i])));
      max_scale= std::max(max_scale, static_cast<Scalar>(atol_ + fabs(scale[i]) * rtol_));
    }
    return std::make_tuple(max_diff, max_scale);
  }
}
#endif //SPACEHUB_IAS15_ERROR_HPP
