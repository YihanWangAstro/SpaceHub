//
// Created by 王艺涵 on 10/9/19.
//

#ifndef SPACEHUB_PI_CONTROLLER_H
#define SPACEHUB_PI_CONTROLLER_H

#include "../stepsize-controller.h"
#include "../own-math.hpp"

namespace space {

  template<size_t Max_order, typename T>
  class PIController : public StepController<PIController<Max_order, T>> {
  public:
    //Type member
    using Base = StepController<PIController<Max_order, T>>;

    using Scalar = T;

    using value_type = T;

    // Constructors
    //SPACEHUB_MAKE_CONSTRUCTORS(PIController, delete, default, default, default, default);

    explicit PIController() {
      for (size_t i = 1; i <= Max_order; i++) {
        expon_[i] = 1.0 / static_cast<Scalar>(i);
        step_limiter_max_[i] = pow(1.0 / safe_factor3_, expon_[i]);
        step_limiter_min_[i] = pow(safe_factor3_, expon_[i]) / safe_factor4_;
      }
    }

    CRTP_impl :
    // CRTP implementation

    template <typename ArrayLike>
    Scalar impl_next_step_size(size_t order, Scalar old_step, ArrayLike const& errors);

  private:
    constexpr static size_t max_order{Max_order};

    std::array<Scalar, Max_order + 1> step_limiter_max_;

    std::array<Scalar, Max_order + 1> step_limiter_min_;

    std::array<Scalar, Max_order + 1> expon_;

    Scalar safe_factor1_{0.94};

    Scalar safe_factor2_{0.65};

    Scalar safe_factor3_{0.02};

    Scalar safe_factor4_{4.0};

    Scalar alpha_{0.7};

    Scalar beta_{0.4};//Proportional Integration feedback
  };

  template<size_t Max_order, typename T>
  template <typename ArrayLike>
  auto PIController<Max_order, T>::impl_next_step_size(size_t order, Scalar old_step, ArrayLike const& errors) -> Scalar {
    static_assert(std::tuple_size<ArrayLike>::value >= 2, "At least 2 numbers of errors need for Proportional Integral controller!");
    if (std::get<0>(errors) != 0.0) {
      return old_step * space::in_range(step_limiter_min_[order],
                                        safe_factor1_ * pow(safe_factor2_ / std::get<0>(errors), alpha_ * expon_[order]) * pow(std::get<1>(errors), beta_ * expon_[order]),
                                        step_limiter_max_[order]);
    } else {
      return old_step * step_limiter_max_[order];
    }
  }
}
#endif //SPACEHUB_PI_CONTROLLER_H
