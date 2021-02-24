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
    the terms of the GPL-3.0 License. SpaceHub is distributed in the hope that it
    will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GPL-3.0 License
    for more details. You should have received a copy of the GPL-3.0 License along
    with SpaceHub.
\*---------------------------------------------------------------------------*/
/**
 * @file Bulirsch-Stoer.hpp
 *
 * Header file.
 */
#pragma once

#include <array>
#include <functional>
#include <vector>

#include "../core-computation.hpp"
#include "../integrator/symplectic/symplectic-integrator.hpp"
#include "../spacehub-concepts.hpp"

#ifdef BS_SUBSTEP_PARALLEL
#include "../taskflow/taskflow.hpp"
#endif

namespace space::ode_iterator {

    /*---------------------------------------------------------------------------*\
         Class BulirschStoerConsts Declaration
    \*---------------------------------------------------------------------------*/
    /**
     * @brief
     *
     * @tparam T
     * @tparam MaxIter
     * @tparam IsKDK
     */
    template <typename T, size_t MaxIter, bool IsKDK>
    class BulirschStoerConsts {
       public:
        using Scalar = T;

        constexpr static size_t max_iter{MaxIter};

        constexpr static double dec_factor{0.8};

        constexpr static double inc_factor{0.9};

        constexpr inline Scalar cost(size_t i) const { return cost_[i]; };

        constexpr inline size_t h(size_t i) const { return h_[i]; };

        constexpr inline Scalar table_coef(size_t i, size_t j) const { return extrap_coef_[at(i, j)]; };

        constexpr explicit BulirschStoerConsts();

       private:
        inline size_t at(size_t i, size_t j) const { return i * MaxIter + j; };

        /** @brief Extrapolation coefficient.*/
        std::array<Scalar, MaxIter *(MaxIter)> extrap_coef_;

        /** @brief The work(computation resource) per step size of each iteration depth.*/
        std::array<Scalar, MaxIter> cost_;

        /** @brief Steps of integration of each iteration depth.*/
        std::array<size_t, MaxIter> h_;
    };

    /*---------------------------------------------------------------------------*\
         Class BulirschStoer Declaration
    \*---------------------------------------------------------------------------*/
    /**
     * @brief
     *
     * @tparam Integrator
     * @tparam ErrEstimator
     * @tparam StepController
     */
    template <typename Integrator, typename ErrEstimator, typename StepController, size_t MaxIter = 7>
    class BulirschStoer {
       public:
        static constexpr size_t max_depth{MaxIter};

        SPACEHUB_USING_TYPE_SYSTEM_OF(Integrator);

        using BSConsts =
            BulirschStoerConsts<Scalar, max_depth + 1, std::is_same_v<Integrator, integrator::LeapFrogKDK<TypeSet>>>;

        using EvaluateFun = std::function<void(StateScalarArray const &, StateScalarArray &, Scalar)>;

        static constexpr size_t max_try_num{100};

        BulirschStoer();

        Scalar iterate(concepts::ParticleSystem auto &particles, Scalar macro_step_size);

        Scalar iterate(EvaluateFun func, StateScalarArray &data, Scalar &time, Scalar step_size);

        void set_atol(Scalar atol);

        void set_rtol(Scalar rtol);

        Scalar reject_rate() { return static_cast<Scalar>(rej_num_) / static_cast<Scalar>(iter_num_); };

       private:
        void check_variable_size();

        inline size_t at(size_t i, size_t j) const { return consts_.at(i, j); };

        void integrate_by_n_steps(concepts::ParticleSystem auto &particles, Scalar macro_step_size, size_t steps);

        void integrate_by_n_steps(EvaluateFun func, StateScalarArray &data_out, Scalar &time, Scalar step_size,
                                  size_t steps);

        void extrapolate(size_t k);

        inline bool in_converged_window(size_t k);

        inline size_t allowed(size_t i) const;

        Scalar set_next_iteration(size_t k);

        bool is_diverged_anyhow(Scalar error, size_t k) const;

        Scalar get_next_step_len(size_t k_new, size_t k) const;

       private:
        /** @brief The constant coef for BS extrapolation*/
        BSConsts consts_;

        /** @brief Extrapolation table.*/
        std::array<ScalarArray, max_depth + 1> extrap_list_;

        ErrEstimator err_checker_;

        StepController step_ctrl_;

        StateScalarArray input_{0};

        /** @brief The optimal step size array.*/
        std::array<Scalar, max_depth + 1> ideal_step_size_{0};

        /** @brief The work(computation resource) needed to converge at column k.*/
        std::array<Scalar, max_depth + 1> cost_per_len_{0};

        Scalar last_error_{1.0};

        /** @brief Current iteration depth.*/
        size_t ideal_rank_{MaxIter - 1};  // make sure the first step reach maxiter

        /** @brief Total volume of extrapolation table(in scalar).*/
        size_t var_num_{0};

        /** @brief Rejection number*/
        size_t rej_num_{0};

        /** @brief Total iteration number*/
        size_t iter_num_{0};

        bool step_reject_{false};

        bool first_step_{true};

#ifdef BS_SUBSTEP_PARALLEL
        tf::Executor executor{1};
#endif
    };

#define CLASS_BulirschStoerConsts(...)                     \
    template <typename Scalar, size_t MaxIter, bool IsKDK> \
    __VA_ARGS__ BulirschStoerConsts<Scalar, MaxIter, IsKDK>

    /*---------------------------------------------------------------------------*\
         Class BulirschStoerConsts Implementation
    \*---------------------------------------------------------------------------*/
    CLASS_BulirschStoerConsts(constexpr)::BulirschStoerConsts() {
        std::array<size_t, 11> seq = {1, 2, 3, 5, 8, 12, 17, 25, 36, 51, 73};  // better sequence

        for (size_t i = 0; i < MaxIter; ++i) {
            // h_[i] = 2 * (i + 1);

            if constexpr (MaxIter <= 11) {
                h_[i] = seq[i];
            } else {
                h_[i] = (i + 1);
            }
            if (i == 0) {
                cost_[i] = h_[i] + static_cast<size_t>(IsKDK) * 2;
            } else {
                cost_[i] = cost_[i - 1] + h_[i] + static_cast<size_t>(IsKDK) * 2;
            }

            for (size_t j = 0; j < MaxIter; ++j) {
                if (j < i) {
                    Scalar nj2 = static_cast<Scalar>(h_[i - j - 1] * h_[i - j - 1]);
                    Scalar ni2 = static_cast<Scalar>(h_[i] * h_[i]);
                    extrap_coef_[at(i, j)] = nj2 / (ni2 - nj2);
                } else {
                    extrap_coef_[at(i, j)] = 0;
                }
            }
        }
    }

#define CLASS_BulirschStoer(...)                                                                   \
    template <typename Integrator, typename ErrEstimator, typename StepController, size_t MaxIter> \
    __VA_ARGS__ BulirschStoer<Integrator, ErrEstimator, StepController, MaxIter>

    /*---------------------------------------------------------------------------*\
         Class BulirschStoer Implementation
    \*---------------------------------------------------------------------------*/
    CLASS_BulirschStoer()::BulirschStoer() : step_ctrl_{}, err_checker_{0, 1e-14} {
        step_ctrl_.set_limiter(0.02, 4.0);
        if (MaxIter < 11) {
            step_ctrl_.set_safe_guards(0.72, 0.95);  // for standard double precision
        } else {
            step_ctrl_.set_safe_guards(0.6, 0.95);  // for arbitrary bits floating points
        }
    }

    CLASS_BulirschStoer(auto)::iterate(EvaluateFun func, StateScalarArray &data, Scalar &time, Scalar step_size)
        -> Scalar {
        Scalar iter_h = step_size;
        input_ = data;

        check_variable_size();

        for (size_t i = 0; i < max_try_num; ++i) {
            iter_num_++;

            integrate_by_n_steps(func, extrap_list_[0], time, iter_h, consts_.h(0));

            for (size_t k = 1; k <= ideal_rank_ + 1; ++k) {
                size_t result_order = 2 * k + 1;
                integrate_by_n_steps(func, extrap_list_[k], time, iter_h, consts_.h(k));

                extrapolate(k);

                Scalar error = err_checker_.error(input_, extrap_list_[1], extrap_list_[0]);

                ideal_step_size_[k] = iter_h * step_ctrl_.next(result_order, error);

                cost_per_len_[k] = consts_.cost(k) / ideal_step_size_[k];
                // space::print_csv(std::cout, k, ideal_rank_, error, ideal_step_size_[k], cost_per_len_[k], '\n');
                if (in_converged_window(k)) {
                    if (error <= 1.0) {
                        step_reject_ = false;
                        time += iter_h;
                        data = extrap_list_[0];
                        Scalar new_h = set_next_iteration(k);
                        iter_h *= step_ctrl_.limiter(result_order, new_h / iter_h);
                        last_error_ = error;
                        return iter_h;
                    } else if (is_diverged_anyhow(error, k)) {
                        step_reject_ = true;
                        rej_num_++;
                        Scalar new_h = set_next_iteration(k);
                        iter_h *= step_ctrl_.limiter(result_order, new_h / iter_h);
                        break;
                    }
                }
            }
            // particles.read_from_scalar_array(input_);
        }
        spacehub_abort("Reach max iteration loop number!");
    }

    CLASS_BulirschStoer(auto)::iterate(concepts::ParticleSystem auto &particles, Scalar macro_step_size) -> Scalar {
        Scalar iter_h = macro_step_size;
        particles.write_to_scalar_array(input_);
        check_variable_size();
        auto const &dy = particles.increment();
        particles.collect_increment(true);

#ifdef BS_SUBSTEP_PARALLEL
        static const auto integrate_with_stride = [&](auto &ptc, size_t start, size_t stride) {
            auto const &dy = ptc.increment();

            for (int p = start; p >= 0; p -= stride) {
                ptc.read_from_scalar_array(input_);
                ptc.clear_increment();
                integrate_by_n_steps(ptc, iter_h, consts_.h(p));
                std::copy(dy.begin(), dy.end(), extrap_list_[p].begin());
            }
        };
#endif

        for (size_t i = 0; i < max_try_num; ++i) {
            iter_num_++;

#ifdef BS_SUBSTEP_PARALLEL
            executor.silent_async(integrate_with_stride, particles, ideal_rank_, 2);
            integrate_with_stride(particles, ideal_rank_ + 1, 2);
            executor.wait_for_all();
#else
            particles.clear_increment();
            integrate_by_n_steps(particles, iter_h, consts_.h(0));
            std::copy(dy.begin(), dy.end(), extrap_list_[0].begin());
#endif

            for (size_t k = 1; k <= ideal_rank_ + 1; ++k) {
                size_t result_order = 2 * k + 1;
#ifdef BS_SUBSTEP_PARALLEL

#else
                particles.read_from_scalar_array(input_);
                particles.clear_increment();
                integrate_by_n_steps(particles, iter_h, consts_.h(k));
                std::copy(dy.begin(), dy.end(), extrap_list_[k].begin());
#endif
                extrapolate(k);  // extrapolate results and save it to extrap_list_[0];
                Scalar error = err_checker_.error(input_, extrap_list_[1], extrap_list_[0]);
                ideal_step_size_[k] = iter_h * step_ctrl_.next(result_order, error);
                cost_per_len_[k] = consts_.cost(k) / ideal_step_size_[k];
                // space::print_csv(std::cout, k, ideal_rank_, error, ideal_step_size_[k], cost_per_len_[k], '\n');
                if (in_converged_window(k)) {
                    if (error <= 1.0) {
                        last_error_ = error;
                        step_reject_ = false;
                        calc::array_advance(input_, extrap_list_[0]);
                        particles.read_from_scalar_array(input_);
                        particles.collect_increment(false);  // turn off the increment collection

                        Scalar new_h = set_next_iteration(k);
                        first_step_ = false;
                        return iter_h * step_ctrl_.limiter(result_order, new_h / iter_h);
                    } else if (is_diverged_anyhow(error, k)) {
                        step_reject_ = true;
                        rej_num_++;
                        Scalar new_h = set_next_iteration(k);
                        iter_h *= step_ctrl_.limiter(result_order, new_h / iter_h);
                        break;
                    }
                }
            }
            particles.read_from_scalar_array(input_);
        }
        spacehub_abort("Reach max iteration loop number!");
    }

    CLASS_BulirschStoer(void)::set_atol(Scalar atol) { err_checker_.set_atol(atol); }

    CLASS_BulirschStoer(void)::set_rtol(Scalar rtol) { err_checker_.set_rtol(rtol); }

    CLASS_BulirschStoer(void)::check_variable_size() {
        var_num_ = input_.size();
        if (var_num_ > extrap_list_[0].size()) [[unlikely]] {
            for (auto &v : extrap_list_) {
                v.resize(var_num_);
                calc::array_set_zero(v);
            }
        }
    }

    CLASS_BulirschStoer(void)::integrate_by_n_steps(concepts::ParticleSystem auto &particles, Scalar macro_step_size,
                                                    size_t steps) {
        Scalar h = macro_step_size / steps;

        if constexpr (std::is_same_v<Integrator, integrator::LeapFrogDKD<TypeSet>>) {
            particles.drift(0.5 * h);
            for (size_t i = 1; i < steps; i++) {
                particles.kick(h);
                particles.drift(h);
            }
            particles.kick(h);
            particles.drift(0.5 * h);
        } else if constexpr (std::is_same_v<Integrator, integrator::LeapFrogKDK<TypeSet>>) {
            particles.kick(0.5 * h);
            for (size_t i = 1; i < steps; i++) {
                particles.drift(h);
                particles.kick(h);
            }
            particles.drift(h);
            particles.kick(0.5 * h);
        } else {
            static_assert(true, "Bulirsch-Stoer Undefined embedded integration method");
        }
    }

    CLASS_BulirschStoer(void)::integrate_by_n_steps(EvaluateFun func, StateScalarArray &data_out, Scalar &time,
                                                    Scalar step_size, size_t steps) {
        data_out = input_;
        Scalar h = step_size / steps;
        StateScalarArray dxdt(input_.size());
        StateScalarArray data_mid(input_);

        func(input_, dxdt, time);
        calc::array_advance(data_out, input_, dxdt, h);
        time += h;
        for (size_t i = 1; i < steps; i++) {
            func(data_out, dxdt, time);
            calc::array_advance(data_mid, dxdt, 2 * h);
            time += h;
            std::swap(data_mid, data_out);
        }
        func(data_out, dxdt, time);
        calc::array_advance(data_mid, dxdt, h);
        calc::array_add(data_out, data_mid, data_out);
        calc::array_scale(data_out, data_out, 0.5);
    }

    CLASS_BulirschStoer(void)::extrapolate(size_t k) {
        for (size_t j = k; j > 0; --j) {
            for (size_t i = 0; i < var_num_; ++i) {
                extrap_list_[j - 1][i] =
                    extrap_list_[j][i] + (extrap_list_[j][i] - extrap_list_[j - 1][i]) * consts_.table_coef(k, k - j);
            }
        }
    }

    CLASS_BulirschStoer(bool)::in_converged_window(size_t k) {
        return (k == ideal_rank_ - 1 || k == ideal_rank_ || k == ideal_rank_ + 1) || (first_step_);
    }

    CLASS_BulirschStoer(size_t)::allowed(size_t i) const {
        return math::in_range(static_cast<size_t>(2), i, static_cast<size_t>(max_depth - 1));
    }

    CLASS_BulirschStoer(auto)::get_next_step_len(size_t k_new, size_t k) const -> Scalar {
        if (k_new <= k) {
            return ideal_step_size_[k_new];
        } else if (k_new == k + 1) {
            return ideal_step_size_[k] * static_cast<Scalar>(consts_.cost(k + 1)) /
                   static_cast<Scalar>(consts_.cost(k));
        } else {
            spacehub_abort("unexpected order!");
        }
    }

    CLASS_BulirschStoer(auto)::set_next_iteration(size_t k) -> Scalar {
        if (!first_step_) {
            if (k == ideal_rank_ - 1 || k == ideal_rank_) [[likely]] {
                if (cost_per_len_[k - 1] < BSConsts::dec_factor * cost_per_len_[k]) [[unlikely]] {
                    ideal_rank_ = allowed(k - 1);
                    return get_next_step_len(ideal_rank_, k);
                } else if (cost_per_len_[k] < BSConsts::inc_factor * cost_per_len_[k - 1] && !step_reject_)
                    [[unlikely]] {
                    ideal_rank_ = allowed(k + 1);
                    return get_next_step_len(ideal_rank_, k);
                } else {
                    ideal_rank_ = allowed(k);
                    return get_next_step_len(ideal_rank_, k);
                }
            } else if (k == ideal_rank_ + 1) {
                if (cost_per_len_[k - 2] < BSConsts::dec_factor * cost_per_len_[k - 1]) {
                    ideal_rank_ = allowed(k - 2);
                }
                if (cost_per_len_[k] < BSConsts::inc_factor * cost_per_len_[ideal_rank_] && !step_reject_) {
                    ideal_rank_ = allowed(k);
                }
                return get_next_step_len(ideal_rank_, k);
            } else {
                spacehub_abort("unexpected iteration index!");
            }
        } else {
            if (!step_reject_) {
                ideal_rank_ = k;
            }
            return ideal_step_size_[k];
        }
    }

    CLASS_BulirschStoer(bool)::is_diverged_anyhow(Scalar error, size_t k) const {
        if (!first_step_) {
            Scalar r = 1.0;
            if (k == ideal_rank_ - 1) {
                r = static_cast<Scalar>(consts_.h(k + 1) * consts_.h(k + 2)) /
                    static_cast<Scalar>(consts_.h(0) * consts_.h(0));
            } else if (k == ideal_rank_) {
                r = static_cast<Scalar>(consts_.h(k + 1)) / static_cast<Scalar>(consts_.h(0));
            }  // else k == iterDepth+1 and error >1 reject directly
            return error > r * r;
        } else {
            return false;
        }
    }

}  // namespace space::ode_iterator
