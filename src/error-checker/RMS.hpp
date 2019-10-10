//
// Created by root on 10/8/19.
//

#ifndef SPACEHUB_RMS_HPP
#define SPACEHUB_RMS_HPP

#include "../error-checker.h"
#include "../own-math.hpp"

namespace space {

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

    CRTP_impl :
    // CRTP implementation

    SPACEHUB_STD_ACCESSOR(auto, impl_atol, atol_);

    SPACEHUB_STD_ACCESSOR(auto, impl_rtol, rtol_);

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
  void RMS<T>::impl_set_atol(Scalar error) {
    atol_ = error;
  }

  template<typename T>
  void RMS<T>::impl_set_rtol(Scalar error) {
    rtol_ = error;
  }

  template<typename T>
  template<typename Array>
  auto RMS<T>::impl_error(const Array &y0, const Array &y1) -> typename Array::value_type {
    size_t const size = y0.size();
    Scalar error = 0;
    for (size_t i = 0; i < size; ++i) {
      auto r = fabs(y0[i] - y1[i]) / (atol_ + std::max(fabs(y0[i]), fabs(y1[i])) * rtol_);
      error += r * r;
    }
    return sqrt(error / size);
  }

  template<typename T>
  template<typename Array>
  auto RMS<T>::impl_error(const Array &scale, const Array &y0, const Array &y1) -> typename Array::value_type {
    size_t const size = scale.size();
    Scalar error = 0;
    for (size_t i = 0; i < size; ++i) {
      auto r = fabs(y0[i] - y1[i]) / (atol_ + fabs(scale[i]) * rtol_);
      error += r * r;
    }
    return sqrt(error / size);
  }
}
#endif //SPACEHUB_RMS_HPP
