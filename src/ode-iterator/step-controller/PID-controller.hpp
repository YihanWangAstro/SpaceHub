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
 * @file PID-controller.hpp
 *
 * Header file.
 */
#ifndef SPACEHUB_PID_CONTROLLER_HPP
#define SPACEHUB_PID_CONTROLLER_HPP

#include "../../math.hpp"
#include "stepsize-controller.hpp"

namespace space::ode_iterator {
/*---------------------------------------------------------------------------*\
     Class PIDController Declaration
\*---------------------------------------------------------------------------*/
/**
 * https://en.wikipedia.org/wiki/PID_controller
 * @tparam Max_order
 * @tparam T
 */
template <size_t Max_order, typename T>
class PIDController : public StepController<PIDController<Max_order, T>> {
 public:
  // Type member
  using Base = StepController<PIDController<Max_order, T>>;

  using Scalar = T;

  using value_type = T;

  constexpr static size_t max_order{Max_order};

  // Constructors
  // SPACEHUB_MAKE_CONSTRUCTORS(PIDController, delete, default, default, default, default);

  explicit PIDController();

  void set_PID_coefficients(Scalar Kp, Scalar Ki, Scalar Kd);

  void set_safe_guards(Scalar S1, Scalar S2, Scalar S3, Scalar S4);

  CRTP_IMPL :
      // CRTP implementation

      template <typename ArrayLike>
      Scalar
      impl_next_step_size(size_t order, Scalar old_step, ArrayLike const &errors);

  Scalar impl_next_step_size(size_t order, Scalar old_step, Scalar error);

 private:
  std::array<Scalar, Max_order + 1> limiter_max_;

  std::array<Scalar, Max_order + 1> limiter_min_;

  std::array<Scalar, Max_order + 1> expon_;

  Scalar safe_guard1_{0.94};

  Scalar safe_guard2_{0.65};

  Scalar safe_guard3_{0.02};

  Scalar safe_guard4_{4.0};

  Scalar Kp_{0.7};  // Proportion feedback coefficient

  Scalar Ki_{0.4};  // Integral feedback coefficient

  Scalar Kd_{0};  // Derivative feedback coefficient

  inline Scalar step_limiter(size_t order, Scalar step_size_ratio);
};
/*---------------------------------------------------------------------------*\
     Class PIDController Implementation
\*---------------------------------------------------------------------------*/
template <size_t Max_order, typename T>
inline auto PIDController<Max_order, T>::step_limiter(size_t order, Scalar step_size_ratio) -> Scalar {
  return math::in_range(limiter_min_[order], step_size_ratio, limiter_max_[order]);
}

template <size_t Max_order, typename T>
void PIDController<Max_order, T>::set_PID_coefficients(Scalar Kp, Scalar Ki, Scalar Kd) {
  Kp_ = Kp, Ki_ = Ki, Kd_ = Kd;
}

template <size_t Max_order, typename T>
void PIDController<Max_order, T>::set_safe_guards(Scalar S1, Scalar S2, Scalar S3, Scalar S4) {
  safe_guard1_ = S1;
  safe_guard2_ = S2;
  safe_guard3_ = S3;
  safe_guard4_ = S4;
}

template <size_t Max_order, typename T>
template <typename ArrayLike>
auto PIDController<Max_order, T>::impl_next_step_size(size_t order, Scalar old_step, ArrayLike const &errors)
    -> Scalar {
  if constexpr (std::tuple_size_v<ArrayLike> == 1) {  // Only proportion part is provided
    if (std::get<0>(errors) != 0.0) {
      return old_step * step_limiter(order, safe_guard1_ * pow(safe_guard2_ / std::get<0>(errors), expon_[order]));
    } else {
      return old_step * limiter_max_[order];
    }
  } else if constexpr (std::tuple_size_v<ArrayLike> == 2) {  // Proportion & Integral part are provided
    if (std::get<0>(errors) != 0.0) {
      return old_step *
             step_limiter(order, safe_guard1_ * pow(safe_guard2_ / std::get<0>(errors), Kp_ * expon_[order]) *
                                     pow(std::get<1>(errors), Ki_ * expon_[order]));
    } else {
      return old_step * limiter_max_[order];
    }
  }
}

template <size_t Max_order, typename T>
auto PIDController<Max_order, T>::impl_next_step_size(size_t order, Scalar old_step, Scalar error) -> Scalar {
  if (error != 0.0) {
    return old_step * step_limiter(order, safe_guard1_ * pow(safe_guard2_ / error, expon_[order]));
  } else {
    return old_step * limiter_max_[order];
  }
}

template <size_t Max_order, typename T>
PIDController<Max_order, T>::PIDController() {
  for (size_t i = 1; i <= Max_order; i++) {
    expon_[i] = 1.0 / static_cast<Scalar>(i);
    limiter_max_[i] = pow(1.0 / safe_guard3_, expon_[i]);
    limiter_min_[i] = pow(safe_guard3_, expon_[i]) / safe_guard4_;
  }
}

}  // namespace space::ode_iterator
#endif  // SPACEHUB_PID_CONTROLLER_HPP
