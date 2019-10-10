//
// Created by 王艺涵 on 10/9/19.
//

#ifndef SPACEHUB_LOCAL_CONTROLLER_H
#define SPACEHUB_LOCAL_CONTROLLER_H

#include "../stepsize-controller.h"
#include "../own-math.hpp"

namespace space {

  template<size_t Max_order, typename T>
  class LocalController : public StepController<LocalController<Max_order, T>> {
  public:
    //Type member
    using Base = StepController<LocalController<Max_order, T>>;

    using Scalar = T;

    using value_type = T;

    // Constructors
    SPACEHUB_MAKE_CONSTRUCTORS(LocalController, delete, default, default, default, default);

    explicit LocalController(Scalar S1 = 0.94, Scalar S2 = 0.65, Scalar S3 = 0.02, Scalar S4 = 4.0) : safe_factor1{S1}, safe_factor2{S2}, safe_factor3{S3}, safe_factor4{S4} {
      for(size_t i = 1 ; i <= Max_order; i++){
        expon_[i] = 1.0 / static_cast<Scalar>(i);
        step_limiter_max_[i] = pow(1.0 / safe_factor3, expon_[i]);
        step_limiter_min_[i] = pow(safe_factor3, expon_[i]) / safe_factor4;
      }
    }

    CRTP_impl :
    // CRTP implementation

    Scalar impl_new_step_size(size_t order, Scalar old_step, Scalar error);

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
  auto LocalController<Max_order, T>::impl_new_step_size(size_t order, Scalar old_step, Scalar error) -> Scalar {
    if (error != 0.0) {
      return old_step * space::in_range(step_limiter_min_[order], safe_factor1 * pow(safe_factor2/error, expon_[order]), step_limiter_max_[order]);
    } else {
      return old_step * step_limiter_max_[order];
    }
  }

  }
#endif //SPACEHUB_LOCAL_CONTROLLER_H
