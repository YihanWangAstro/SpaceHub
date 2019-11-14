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
 * @file Burlish-Stoer.hpp
 *
 * Header file.
 */
#ifndef SPACEHUB_BURLISH_STOER_HPP
#define SPACEHUB_BURLISH_STOER_HPP

#include <array>
#include <vector>
#include "ode-iterator.hpp"

namespace space::ode_iterator {
/*---------------------------------------------------------------------------*\
     Class BurlishStoerConsts Declaration
\*---------------------------------------------------------------------------*/
/**
 *
 * @tparam T
 * @tparam MaxIter
 */
template <typename T, size_t MaxIter>
class BurlishStoerConsts {
 public:
  using Scalar = T;

  constexpr static size_t max_iter{MaxIter};

  constexpr static Scalar cost_tol{0.9};

  inline Scalar cost(size_t i) const;

  [[nodiscard]] inline size_t step_sequence(size_t i) const;

  inline Scalar table_coef(size_t i, size_t j) const;

  explicit BurlishStoerConsts();

 private:
  [[nodiscard]] inline size_t at(size_t i, size_t j) const;

  /** @brief Extrapolation coefficient.*/
  std::array<Scalar, MaxIter *(MaxIter)> extrap_coef_;

  /** @brief The work(computation resource) per step size of each iteration depth.*/
  std::array<size_t, MaxIter> cost_;

  /** @brief Steps of integration of each iteration depth.*/
  std::array<size_t, MaxIter> sub_steps_;
};

/*---------------------------------------------------------------------------*\
     Class BurlishStoer Declaration
\*---------------------------------------------------------------------------*/
/**
 *
 * @tparam Real
 * @tparam ErrChecker
 * @tparam StepControl
 */
template <typename Real, template <typename> typename ErrChecker, template <size_t, typename> typename StepControl>
class BurlishStoer : public OdeIterator<BurlishStoer<Real, ErrChecker, StepControl>> {
  static_assert(std::is_floating_point<Real>::value, "Only float-like type can be used!");

 public:
  using Base = OdeIterator<BurlishStoer<Real, ErrChecker, StepControl>>;

  using Scalar = Real;

  static constexpr size_t max_depth{7};

  static constexpr size_t max_try_num{500};

  using BSConsts = BurlishStoerConsts<Scalar, max_depth + 1>;

  using StepController = StepControl<2 * max_depth + 3, Scalar>;

  CRTP_IMPL :

      template <typename U>
      Scalar
      impl_iterate(U &particles, typename U::Scalar macro_step_size);

  void impl_set_atol(Scalar atol);

  void impl_set_rtol(Scalar rtol);

 private:
  void check_variable_size();

  [[nodiscard]] inline size_t at(size_t i, size_t j) const;

  template <typename U>
  void integrate_by_n_steps(U &particles, Scalar macro_step_size, size_t steps);

  void extrapolate(size_t k);

  inline bool in_converged_window(size_t k);

  Scalar get_next_try_step(size_t k);

  [[nodiscard]] inline size_t allowed(size_t i) const;

  Scalar set_next_iteration(size_t k, bool last_reject);

  bool is_diverged_anyhow(Scalar error, size_t k) const;

 private:
  /** @brief The constat coef for BS extrapolation*/
  BSConsts parameters_;

  /** @brief Extrapolation table.*/
  std::array<std::vector<Scalar>, max_depth + 1> extrap_list_;

  ErrChecker<Scalar> err_checker_;

  StepController step_controller_;

  std::vector<Scalar> input_;

  /** @brief The optimal step size array.*/
  std::array<Scalar, max_depth + 1> ideal_step_size_{0};

  /** @brief The work(computation resource) needed to converge at column k.*/
  std::array<Scalar, max_depth + 1> cost_per_len_{0};

  Scalar last_error_{1.0};

  /** @brief Current iteraation depth.*/
  size_t ideal_rank_{4};

  /** @brief Total volume of extrapolation table(in scalar).*/
  size_t var_num_{0};

  /** @brief Rejectin number*/
  size_t rej_num_{0};

  /** @brief Total iteration number*/
  size_t iter_num_{0};

  bool last_step_reject_{false};
};

/*---------------------------------------------------------------------------*\
     Class BurlishStoerConsts Implementation
\*---------------------------------------------------------------------------*/
template <typename Scalar, size_t MaxIter>
Scalar BurlishStoerConsts<Scalar, MaxIter>::cost(size_t i) const {
  return cost_[i];
}

template <typename Scalar, size_t MaxIter>
size_t BurlishStoerConsts<Scalar, MaxIter>::step_sequence(size_t i) const {
  return sub_steps_[i];
}

template <typename Scalar, size_t MaxIter>
Scalar BurlishStoerConsts<Scalar, MaxIter>::table_coef(size_t i, size_t j) const {
  return extrap_coef_[at(i, j)];
}

template <typename Scalar, size_t MaxIter>
BurlishStoerConsts<Scalar, MaxIter>::BurlishStoerConsts() {
  for (size_t i = 0; i < MaxIter; ++i) {
    sub_steps_[i] = 2 * (i + 1);

    if (i == 0) {
      cost_[i] = sub_steps_[i];  // The additional 1 is for 'KDK' method.
    } else {
      cost_[i] = cost_[i - 1] + sub_steps_[i];
    }

    for (size_t j = 0; j < i; ++j) {
      Scalar ratio = static_cast<Scalar>(sub_steps_[i]) / static_cast<Scalar>(sub_steps_[j]);
      extrap_coef_[at(i, j)] = 1.0 / (ratio * ratio - 1.0);
    }
  }
}

template <typename Scalar, size_t MaxIter>
size_t BurlishStoerConsts<Scalar, MaxIter>::at(size_t i, size_t j) const {
  // return i * (i + 1) / 2 + j;
  return i * MaxIter + j;
}

/*---------------------------------------------------------------------------*\
     Class BurlishStoer Implementation
\*---------------------------------------------------------------------------*/
template <typename Real, template <typename> typename ErrChecker, template <size_t, typename> typename StepControl>
template <typename U>
auto BurlishStoer<Real, ErrChecker, StepControl>::impl_iterate(U &particles, typename U::Scalar macro_step_size)
    -> Scalar {
  static_assert(particle_system::is_particle_system_v<U>, "Passing non particle-system-type!");

  Scalar iter_h = macro_step_size;
  particles.to_linear_container(input_);
  check_variable_size();

  for (size_t i = 0; i < max_try_num; ++i) {
    iter_num_++;
    bool reject = true;

    integrate_by_n_steps(particles, iter_h, parameters_.step_sequence(0));
    particles.to_linear_container(extrap_list_[0]);

    for (size_t k = 1; k <= ideal_rank_ + 1; ++k) {
      particles.load_from_linear_container(input_);
      integrate_by_n_steps(particles, iter_h, parameters_.step_sequence(k));
      particles.to_linear_container(extrap_list_[k]);
      extrapolate(k);

      Scalar error = err_checker_.error(input_, extrap_list_[0], extrap_list_[1]);

      ideal_step_size_[k] = step_controller_.next_step_size(2 * k + 1, iter_h, error);

      cost_per_len_[k] = parameters_.cost(k) / ideal_step_size_[k];
      // space::print_csv(std::cout, k, ideal_rank_, error, ideal_step_size_[k], cost_per_len_[k],'\n');
      if (in_converged_window(k)) {
        if (error < 1.0) {
          reject = false;
          iter_h = set_next_iteration(k, last_step_reject_);
          particles.load_from_linear_container(extrap_list_[0]);
          last_step_reject_ = reject;
          last_error_ = error;
          return iter_h;
        } else if (is_diverged_anyhow(error, k)) {
          reject = true;
          rej_num_++;
          iter_h = get_next_try_step(k);
          break;
        }
      }
    }
    last_step_reject_ = reject;
    particles.load_from_linear_container(input_);
  }

  spacehub_abort("Reach max iteration loop number!");
}

template <typename Real, template <typename> typename ErrChecker, template <size_t, typename> typename StepControl>
void BurlishStoer<Real, ErrChecker, StepControl>::impl_set_atol(Scalar atol) {
  err_checker_.set_atol(atol);
}

template <typename Real, template <typename> typename ErrChecker, template <size_t, typename> typename StepControl>
void BurlishStoer<Real, ErrChecker, StepControl>::impl_set_rtol(Scalar rtol) {
  err_checker_.set_rtol(rtol);
}

template <typename Real, template <typename> typename ErrChecker, template <size_t, typename> typename StepControl>
void BurlishStoer<Real, ErrChecker, StepControl>::check_variable_size() {
  var_num_ = input_.size();
  if (var_num_ > extrap_list_[0].size()) {
    for (auto &v : extrap_list_) {
      v.resize(var_num_);
    }
  }
}

template <typename Real, template <typename> typename ErrChecker, template <size_t, typename> typename StepControl>
size_t BurlishStoer<Real, ErrChecker, StepControl>::at(size_t i, size_t j) const {
  // return i * (i + 1) / 2 + j;
  return i * BSConsts::max_iter + j;
}

template <typename Real, template <typename> typename ErrChecker, template <size_t, typename> typename StepControl>
template <typename U>
void BurlishStoer<Real, ErrChecker, StepControl>::integrate_by_n_steps(U &particles, Scalar macro_step_size,
                                                                       size_t steps) {
  Scalar h = macro_step_size / steps;
  particles.drift(0.5 * h);
  for (size_t i = 1; i < steps; i++) {
    particles.kick(h);
    particles.drift(h);
  }
  particles.kick(h);
  particles.drift(0.5 * h);
}

template <typename Real, template <typename> typename ErrChecker, template <size_t, typename> typename StepControl>
void BurlishStoer<Real, ErrChecker, StepControl>::extrapolate(size_t k) {
  for (size_t j = k; j > 0; --j) {
    //auto c_1 = 1 + parameters_.table_coef(k, j - 1);
    //auto c_2 = -parameters_.table_coef(k, j - 1);
    for (size_t i = 0; i < var_num_; ++i) {
      //extrap_list_[j - 1][i] = c_1 * extrap_list_[j][i] + c_2 * extrap_list_[j - 1][i];
        extrap_list_[j - 1][i] = extrap_list_[j][i] + (extrap_list_[j][i] - extrap_list_[j - 1][i]) * parameters_.table_coef(k, j - 1);
    }
  }
}

template <typename Real, template <typename> typename ErrChecker, template <size_t, typename> typename StepControl>
bool BurlishStoer<Real, ErrChecker, StepControl>::in_converged_window(size_t k) {
  return k == ideal_rank_ - 1 || k == ideal_rank_ || k == ideal_rank_ + 1;
}

template <typename Real, template <typename> typename ErrChecker, template <size_t, typename> typename StepControl>
auto BurlishStoer<Real, ErrChecker, StepControl>::get_next_try_step(size_t k) -> Scalar {
  if (k == ideal_rank_ - 1 || k == ideal_rank_) {
    return ideal_step_size_[k];
  } else if (k == ideal_rank_ + 1) {
    return ideal_step_size_[ideal_rank_];
  } else {
    spacehub_abort("unexpected iteration index!");
  }
}

template <typename Real, template <typename> typename ErrChecker, template <size_t, typename> typename StepControl>
size_t BurlishStoer<Real, ErrChecker, StepControl>::allowed(size_t i) const {
  return math::in_range(static_cast<size_t>(2), i, static_cast<size_t>(max_depth - 1));
}

template <typename Real, template <typename> typename ErrChecker, template <size_t, typename> typename StepControl>
auto BurlishStoer<Real, ErrChecker, StepControl>::set_next_iteration(size_t k, bool last_reject) -> Scalar {
  if (k == ideal_rank_) {
    if (cost_per_len_[k - 1] < BSConsts::cost_tol * cost_per_len_[k]) {
      ideal_rank_ = allowed(k - 1);
      return ideal_step_size_[ideal_rank_];
    } else if (cost_per_len_[k] < BSConsts::cost_tol * cost_per_len_[k - 1] && !last_reject) {
      ideal_rank_ = allowed(k + 1);
      return ideal_step_size_[k] * static_cast<Scalar>(parameters_.cost(ideal_rank_)) /
             static_cast<Scalar>(parameters_.cost(k));
    } else {
      return ideal_step_size_[ideal_rank_];
    }
  } else if (k == ideal_rank_ - 1) {
    if (ideal_rank_ <= 2 || cost_per_len_[k] < BSConsts::cost_tol * cost_per_len_[k - 1]) {
      ideal_rank_ = allowed(k + 1);
      return ideal_step_size_[k] * static_cast<Scalar>(parameters_.cost(k + 1)) /
             static_cast<Scalar>(parameters_.cost(k));
    } else {
      ideal_rank_ = allowed(k);
      return ideal_step_size_[k];
    }
  } else if (k == ideal_rank_ + 1) {
    if (cost_per_len_[k - 2] < BSConsts::cost_tol * cost_per_len_[k - 1]) {
      ideal_rank_ = allowed(ideal_rank_ - 1);
    }
    if (cost_per_len_[k] < BSConsts::cost_tol * cost_per_len_[ideal_rank_] && !last_reject) {
      ideal_rank_ = allowed(k);
    }
    return ideal_step_size_[ideal_rank_];
  } else {
    spacehub_abort("unexpected iteration index!");
  }
}

template <typename Real, template <typename> typename ErrChecker, template <size_t, typename> typename StepControl>
bool BurlishStoer<Real, ErrChecker, StepControl>::is_diverged_anyhow(Scalar error, size_t k) const {
  Scalar r = 1.0;

  if (k == ideal_rank_ - 1) {
    r = static_cast<Scalar>(parameters_.step_sequence(k + 1) * parameters_.step_sequence(k + 2)) /
        static_cast<Scalar>(parameters_.step_sequence(0) * parameters_.step_sequence(0));
  } else if (k == ideal_rank_) {
    r = static_cast<Scalar>(parameters_.step_sequence(k)) / static_cast<Scalar>(parameters_.step_sequence(0));
  }  // else k == iterDepth+1 and error >1 reject directly

  return error > r * r;
}
}  // namespace space::ode_iterator
#endif  // SPACEHUB_BURLISH_STOER_HPP
