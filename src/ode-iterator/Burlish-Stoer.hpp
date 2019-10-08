//
// Created by 王艺涵 on 10/7/19.
//

#ifndef SPACEHUB_BS_ITERATOR_NEW_H
#define SPACEHUB_BS_ITERATOR_NEW_H

#include <array>
#include <vector>
#include "ode-iterator.hpp"
#include "../own-math.hpp"
#include "../core-computation.hpp"

namespace space::odeIterator {

  template<typename T, size_t MaxIter>
  class BurlishStoerConsts {
  public:
    using Scalar = T;

    constexpr static size_t max_iter{MaxIter};

    constexpr static Scalar S1{0.65};

    constexpr static Scalar S2{0.94};

    constexpr static Scalar S3{0.02};

    constexpr static Scalar S4{4.0};

    constexpr static Scalar F1{0.8};

    constexpr static Scalar F2{0.9};

    inline Scalar expon(size_t i) const {
      return err_expon_[i];
    }

    inline Scalar step_limiter(size_t i) const {
      return step_limiter_[i];
    }

    inline Scalar work(size_t i) const {
      return work_[i];
    }

    inline size_t step_sequence(size_t i) const {
      return sub_steps_[i];
    }

    inline Scalar table_coef(size_t k) const {
      return extrap_coef_[k];
    }

    inline Scalar table_coef(size_t i, size_t j) const {
      return extrap_coef_[at(i, j)];
    }

    explicit BurlishStoerConsts() {
      for (size_t i = 0; i < MaxIter; ++i) {
        sub_steps_[i] = 2 * (i + 1);

        if (i == 0) {
          work_[i] = sub_steps_[i] + 1;//The additional 1 is for 'KDK' method.
        } else {
          work_[i] = work_[i - 1] + sub_steps_[i] + 1;
        }

        err_expon_[i] = 1.0 / (static_cast<Scalar>(i) * 2 + 1);
        step_limiter_[i] = pow(S3, err_expon_[i]);

        for (size_t j = 0; j < i; ++j) {
          Scalar ratio = static_cast<Scalar>(sub_steps_[i]) / static_cast<Scalar>(sub_steps_[j]);
          extrap_coef_[at(i, j)] = 1.0 / (ratio * ratio - 1.0);
        }
      }
    }

    friend std::ostream &operator<<(std::ostream &os, BurlishStoerConsts const &tab) {
      os << tab.sub_steps_ << "\n\n";
      os << tab.work_ << "\n\n";
      os << tab.err_expon_ << "\n\n";
      os << tab.step_limiter_ << "\n\n";
      for (size_t i = 0; i < BurlishStoerConsts::max_iter; ++i) {
        for (size_t j = 0; j < i; ++j) {
          os << tab.extrap_coef_[tab.at(i, j)] << ',';
        }
        os << '\n';
      }
      return os;
    }

  private:
    inline size_t at(size_t i, size_t j) const {
      return i * (i + 1) / 2 + j;
    }

    /** @brief Extrapolation coefficient.*/
    std::array<Scalar, MaxIter * (MaxIter + 1) / 2> extrap_coef_;

    /** @brief The exponent of error estimate at column k.*/
    std::array<Scalar, MaxIter> err_expon_;

    /** @brief The minimal coeffient of integration step estimation.*/
    std::array<Scalar, MaxIter> step_limiter_;

    /** @brief The work(computation resource) per step size of each iteration depth.*/
    std::array<size_t, MaxIter> work_;

    /** @brief Steps of integration of each iteration depth.*/
    std::array<size_t, MaxIter> sub_steps_;
  };

  template<typename Real, template<typename> typename ErrChecker>
  class BurlishStoer : public OdeIterator<BurlishStoer<Real, ErrChecker>> {
    static_assert(std::is_floating_point<Real>::value, "Only float-like type can be used!");
  public:
    using Scalar = Real;

    static constexpr size_t max_depth{8};

    static constexpr size_t max_try_num{500};

    using BSConsts = BurlishStoerConsts<Scalar, max_depth + 1>;

    template<typename U>
    auto impl_iterate(U &ptcs, typename U::Scalar macro_step_size) -> typename U::Scalar {
      static_assert(is_particle_system_v<U>, "Passing non paritcle-system-type!");

      Scalar iter_H = macro_step_size;
      ptcs.to_linear_container(input_);
      check_variable_size();

      for (size_t i = 0; i < max_try_num; ++i) {
        iter_num_++;
        bool reject = true;

        evolve_by_n_steps(ptcs, iter_H, BS_.step_sequence(0));
        ptcs.to_linear_container(extrap_tab_[0]);

        for (size_t k = 1; k <= ideal_iter_ + 1; ++k) {
          ptcs.load_from_linear_container(input_);
          evolve_by_n_steps(ptcs, iter_H, BS_.step_sequence(k));
          ptcs.to_linear_container(extrap_tab_[k]);
          extrapolate_tab(k);
          calc::array_sub(diff_, extrap_tab_[0], extrap_tab_[1]);
          Scalar error = err_checker_.error(input_, diff_);

          optimal_step_size_[k] = calc_ideal_step_coef(iter_H, error, k);

          work_per_len_[k] = BS_.work(k) / optimal_step_size_[k];

          if (in_converged_window(k)) {
            if (error < 1.0) {
              reject = false;
              iter_H = set_next_iteration(k, last_step_reject);
              ptcs.load_from_linear_container(extrap_tab_[0]);
              return iter_H;
            } else if (is_diverged_anyhow(error, k)) {
              reject = true;
              rej_num_++;
              iter_H = get_next_try_step(k);
              break;
            }
          }
        }
        last_step_reject = reject;
        ptcs.load_from_linear_container(input_);
      }

      spacehub_abort("Reach max iteration loop number!");
    }

  private:
    void check_variable_size() {
      var_num_ = input_.size();
      if (var_num_ > diff_.size()) {
        diff_.resize(var_num_);
        for (auto &v : extrap_tab_) {
          v.resize(var_num_);
        }
      }
    }

    inline size_t at(size_t i, size_t j) const {
      return i * (i + 1) / 2 + j;
    }

    template<typename U>
    void evolve_by_n_steps(U &ptcs, Scalar macro_step_size, size_t steps) {
      Scalar h = macro_step_size / steps;
      ptcs.kick(0.5 * h);
      for (size_t i = 1; i < steps; i++) {
        ptcs.drift(h);
        ptcs.kick(h);
      }
      ptcs.drift(h);
      ptcs.kick(0.5 * h);
    }

    void extrapolate_tab(size_t k) {
      for (size_t j = k - 1; j > 0; --j) {
        auto C1 = 1 + BS_.table_coef(k, j);
        auto C2 = -BS_.table_coef(k, j);
        for (size_t i = 0; i < var_num_; ++i) {
          extrap_tab_[j][i] = C1 * extrap_tab_[j + 1][i] + C2 * extrap_tab_[j][i];
        }
      }
      auto C1 = 1 + BS_.table_coef(k, 0);
      auto C2 = -BS_.table_coef(k, 0);
      for (size_t i = 0; i < var_num_; ++i) {
        extrap_tab_[0][i] = C1 * extrap_tab_[1][i] + C2 * extrap_tab_[0][i];
      }
    }

    inline bool in_converged_window(size_t iter) {
      return iter == ideal_iter_ - 1 || iter == ideal_iter_ || iter == ideal_iter_ + 1;
    }

    Scalar get_next_try_step(size_t iter) {
      if (iter == ideal_iter_ - 1 || iter == ideal_iter_) {
        return optimal_step_size_[iter];
      } else if (iter == ideal_iter_ + 1) {
        return optimal_step_size_[ideal_iter_];
      } else {
        spacehub_abort("unexpected iteration index!");
        return optimal_step_size_[iter];
      }
    }

    Scalar calc_ideal_step_coef(Scalar h, Scalar error, size_t k) {
      if (error != 0.0) {
        return h * space::in_range(BS_.step_limiter(k) / BSConsts::S4,
                                   BSConsts::S2 / pow(error / BSConsts::S1, BS_.expon(k)), 1.0 / BS_.step_limiter(k));
      } else {
        return h / BS_.step_limiter(k);
      }
    }

    inline size_t allowed(size_t i) const {
      return space::in_range(static_cast<size_t>(2), i, static_cast<size_t>(max_depth - 1));
    }

    Scalar set_next_iteration(size_t iter, bool last_reject) {
      if (iter == ideal_iter_ - 1) {
        if (ideal_iter_ <= 2 || work_per_len_[iter] < BSConsts::F2 * work_per_len_[iter - 1]) {
          ideal_iter_ = allowed(iter + 1);
          return optimal_step_size_[iter] * static_cast<Scalar>(BS_.work(iter + 1)) /
                 static_cast<Scalar>(BS_.work(iter));
        } else {
          ideal_iter_ = allowed(iter);
          return optimal_step_size_[iter];
        }
      } else if (iter == ideal_iter_) {
        if (work_per_len_[iter - 1] < BSConsts::F2 * work_per_len_[iter]) {
          ideal_iter_ = allowed(iter - 1);
          return optimal_step_size_[ideal_iter_];
        } else if (work_per_len_[iter] < BSConsts::F2 * work_per_len_[iter - 1] && !last_reject) {
          ideal_iter_ = allowed(iter + 1);
          return optimal_step_size_[iter] * static_cast<Scalar>(BS_.work(ideal_iter_)) /
                 static_cast<Scalar>(BS_.work(iter));
        } else {
          return optimal_step_size_[ideal_iter_];
        }
      } else if (iter == ideal_iter_ + 1) {
        if (work_per_len_[iter - 2] < BSConsts::F2 * work_per_len_[iter - 1]) {
          ideal_iter_ = allowed(ideal_iter_ - 1);
        }
        if (work_per_len_[iter] < BSConsts::F2 * work_per_len_[ideal_iter_] && !last_reject) {
          ideal_iter_ = allowed(iter);
        }
        return optimal_step_size_[ideal_iter_];
      } else {
        spacehub_abort("unexpected iteration index!");
      }
    }

    bool is_diverged_anyhow(Scalar error, size_t iter) const {
      Scalar r = 1;

      if (iter == ideal_iter_ - 1) {
        r = static_cast<Scalar>(BS_.step_sequence(iter + 1) * BS_.step_sequence(iter + 2)) /
            static_cast<Scalar>(BS_.step_sequence(0) * BS_.step_sequence(0));
      } else if (iter == ideal_iter_) {
        r = static_cast<Scalar>(BS_.step_sequence(iter)) / static_cast<Scalar>(BS_.step_sequence(0));
      } else {
        return error > 1.0;//k == iterDepth+1 and error >1 reject directly
      }

      return error > r * r;
    }


  private:
    /** @brief The constat coef for BS extrapolation*/
    BSConsts BS_;

    /** @brief Extrapolation table.*/
    std::array<std::vector<Scalar>, max_depth + 1> extrap_tab_;

    ErrChecker<Scalar> err_checker_;

    std::vector<Scalar> diff_;

    std::vector<Scalar> input_;

    /** @brief The optimal step size array.*/
    std::array<Scalar, max_depth + 1> optimal_step_size_{0};

    /** @brief The work(computation resource) needed to converge at column k.*/
    std::array<Scalar, max_depth + 1> work_per_len_{0};

    /** @brief Current iteraation depth.*/
    size_t ideal_iter_{4};

    /** @brief Total volume of extrapolation table(in scalar).*/
    size_t var_num_{0};

    /** @brief Rejectin number*/
    size_t rej_num_{0};

    /** @brief Total iteration number*/
    size_t iter_num_{0};

    bool last_step_reject{false};
  };
}
#endif //SPACEHUB_BS_ITERATOR_NEW_H
