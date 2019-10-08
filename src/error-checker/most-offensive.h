//
// Created by 王艺涵 on 10/8/19.
//

#ifndef SPACEHUB_MOST_OFFENSIVE_H
#define SPACEHUB_MOST_OFFENSIVE_H

#include "../error-checker.h"
#include "../own-math.hpp"
namespace space{

  template <typename T>
  class MostOffensive : public ErrorChecker<MostOffensive<T>> {
  public:
    //Type member
    using Base = ErrorChecker<MostOffensive<T>>;

    using Scalar = T;

    using value_type = T;

    // Constructors
    MostOffensive() = default;

    MostOffensive(Scalar atol, Scalar rtol) : atol_{atol}, rtol_{rtol} {}

    MostOffensive(MostOffensive const &) = default;

    MostOffensive(MostOffensive &&) noexcept = default;

    MostOffensive &operator=(MostOffensive const &) = default;

    MostOffensive &operator=(MostOffensive &&) noexcept = default;

    CRTP_impl :
    // CRTP implementation

    SPACEHUB_STD_ACCESSOR(auto, impl_atol, atol_);

    SPACEHUB_STD_ACCESSOR(auto, impl_rtol, rtol_);

    void impl_set_atol(Scalar);

    void impl_set_rtol(Scalar);

    template <typename Array>
    auto impl_error(Array const& scale, Array const& diff)->typename Array::value_type;

  private:
    Scalar atol_{1e-13};

    Scalar rtol_{1e-13};
  };

  template <typename T>
  void MostOffensive<T>::impl_set_atol(Scalar error) {
    atol_ = error;
  }

  template <typename T>
  void MostOffensive<T>::impl_set_rtol(Scalar error) {
    rtol_ = error;
  }

  template <typename T>
  template <typename Array>
  auto MostOffensive<T>::impl_error(const Array &scale, const Array &diff)->typename Array::value_type {
    size_t const size = scale.size();
    Scalar max_err = 0;
    for(size_t i = 0 ; i < size; ++i){
      max_err = space::max(max_err, fabs(diff[i] / (fabs(scale[i]) + atol_)));
    }
    return max_err/rtol_;
  }
}
#endif //SPACEHUB_MOST_OFFENSIVE_H
