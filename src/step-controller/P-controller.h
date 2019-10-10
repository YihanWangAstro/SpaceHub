//
// Created by 王艺涵 on 10/9/19.
//

#ifndef SPACEHUB_P_CONTROLLER_H
#define SPACEHUB_P_CONTROLLER_H

#include "../stepsize-controller.h"
#include "../own-math.hpp"

namespace space {

  template<size_t Max_order, typename T>
  class PController : public StepController<PController<Max_order, T>> {
  public:
    //Type member
    using Base = StepController<PController<Max_order, T>>;

    using Scalar = T;

    using value_type = T;

    // Constructors
    //SPACEHUB_MAKE_CONSTRUCTORS(PController, default, default, default, default, default);

    explicit PController() {
      for(size_t i = 1 ; i <= Max_order; i++){
        expon_[i] = 1.0 / static_cast<Scalar>(i);
        step_limiter_max_[i] = pow(1.0 / safe_factor3, expon_[i]);
        step_limiter_min_[i] = pow(safe_factor3, expon_[i]) / safe_factor4;
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

    Scalar safe_factor1{0.94};

    Scalar safe_factor2{0.65};

    Scalar safe_factor3{0.02};

    Scalar safe_factor4{4.0};
  };

  template<size_t Max_order, typename T>
  template <typename ArrayLike>
  auto PController<Max_order, T>::impl_next_step_size(size_t order, Scalar old_step, ArrayLike const& errors) -> Scalar {
    if (std::get<0>(errors) != 0.0) {
      return old_step * space::in_range(step_limiter_min_[order], safe_factor1 * pow(safe_factor2/std::get<0>(errors), expon_[order]), step_limiter_max_[order]);
    } else {
      return old_step * step_limiter_max_[order];
    }
  }

  }
#endif //SPACEHUB_P_CONTROLLER_H
