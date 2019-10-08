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
    RMS() = default;

    RMS(Scalar atol, Scalar rtol) : atol_{atol}, rtol_{rtol} {}

    RMS(RMS const &) = default;

    RMS(RMS &&) noexcept = default;

    RMS &operator=(RMS const &) = default;

    RMS &operator=(RMS &&) noexcept = default;

    CRTP_impl :
    // CRTP implementation

    SPACEHUB_STD_ACCESSOR(auto, impl_atol, atol_);

    SPACEHUB_STD_ACCESSOR(auto, impl_rtol, rtol_);

    void impl_set_atol(Scalar);

    void impl_set_rtol(Scalar);

    template<typename Array>
    auto impl_error(Array const &scale, Array const &diff) -> typename Array::value_type;

  private:
    Scalar atol_{1e-13};

    Scalar rtol_{1e-13};
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
    size_t const size = scale.size();
    Scalar error = 0;
    for (size_t i = 0; i < size; ++i) {
      auto r = diff[i] / (atol_ + fabs(scale[i]) * rtol_);
      error += r * r;
    }
    return sqrt(error / size);
  }
}
#endif //SPACEHUB_RMS_HPP
