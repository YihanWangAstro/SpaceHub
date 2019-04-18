//
// Created by yihan on 4/17/19.
//

#ifndef SPACEHUB_BS_ITERATOR0_H
#define SPACEHUB_BS_ITERATOR0_H

#include <array>
#include <vector>
#include "ode-iterator.h"
#include "../own-math.h"

namespace space::odeIterator {
    template<typename T, size_t MaxIter>
    class BStab {
    public:
        using Scalar = T;

        inline Scalar expon(size_t i) const {
            return err_expon_[i];
        }

        inline Scalar limiter(size_t i) const {
            return limiter_[i];
        }

        inline Scalar work(size_t i) const {
            return work_[i];
        }

        inline size_t step(size_t i) const {
            return sub_steps_[i];
        }

        inline Scalar coef(size_t k) const {
            return extrap_coef_[k];
        }

        inline Scalar coef(size_t i, size_t j) const {
            return extrap_coef_[at(i, j)];
        }

        BStab() {
            for (size_t i = 0; i < MaxIter; ++i) {
                if (i == 0) {
                    sub_steps_[i] = 1;
                    work_[i] = 2;
                } else {
                    sub_steps_[i] = 2 * i;
                    work_[i] = work_[i - 1] + sub_steps_[i] + 1;
                }
                err_expon_[i] = static_cast<Scalar>(1.0 / (2 * i + 1));
                limiter_[i] = pow(0.02, err_expon_[i]);

                for (size_t j = 0; j < i; ++j) {
                    Scalar ratio = static_cast<Scalar>(sub_steps_[i]) / static_cast<Scalar>(sub_steps_[i - j - 1]);
                    extrap_coef_[at(i, j)] = 1.0 / (ratio * ratio - 1);
                }
            }
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
        std::array<Scalar, MaxIter> limiter_;

        /** @brief The work(computation resource) per step size of each iteration depth.*/
        std::array<size_t, MaxIter> work_;

        /** @brief Steps of integration of each iteration depth.*/
        std::array<size_t, MaxIter> sub_steps_;
    };

    template<typename Real>
    class BSIterator : public OdeIterator<BSIterator<Real>> {
        static_assert(std::is_floating_point<Real>::value, "Only float-like type can be used!");
    public:
        using Scalar = Real;

        template<typename U>
        auto impl_iterate(U &ptcs, typename U::Scalar macro_step_size) -> typename U::Scalar {
            static_assert(is_particle_system_v<U>, "Passing non paritcle-system-type!");

            Scalar iter_H = macro_step_size;

            for (;;) {
                iter_num_++;
                auto local_sys = ptcs;

                evolve_by_n_steps(local_sys, iter_H, BS_.step(0));

                local_sys.to_linear_container(extrap_tab_[at(0, 0)]);

                var_num_ = extrap_tab_[0].size();

                for (size_t i = 1; i <= ideal_iter_ + 1; ++i) {
                    local_sys = ptcs;
                    evolve_by_n_steps(local_sys, iter_H, BS_.step(i));
                    local_sys.to_linear_container(extrap_tab_[at(i, 0)]);

                    extrapolate_tab(i);

                    Scalar error = calc_error_of_row(i);
                    optm_step_coef_[i] = calc_ideal_step_coef(error, i);
                    work_per_len_[i] = BS_.work(i) / optm_step_coef_[i];

                    if (in_converged_window(i)) {
                        if (error < 1.0) {
                            iter_H *= step_coef_limiter(prepare_next_iter(i));
                            ptcs.load_from_linear_container(extrap_tab_[at(i, i)]);
                            return iter_H;
                        } else if (is_diverged_anyhow(error, i)) {
                            rej_num_++;
                            iter_H *= step_coef_limiter(reduced_step_coef(i));
                            break;
                        }
                    }
                }
            }
        }

        void set_abs_err(Scalar err) {
            abs_error_ = err;
        }

        void set_rel_err(Scalar err) {
            rel_error_ = err;
        }

    private:
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
            size_t curr_row = at(k, 0);
            size_t last_row = at(k - 1, 0);
            size_t size = extrap_tab_[curr_row].size();
            for (size_t j = 0; j < k; ++j) {
                size_t center = curr_row + j;//
                size_t right = center + 1;
                size_t up = last_row + j;

                if (extrap_tab_[right].size() != var_num_)
                    extrap_tab_[right].resize(var_num_);

                for (size_t i = 0; i < var_num_; ++i)
                    extrap_tab_[right][i] =
                            extrap_tab_[center][i] + (extrap_tab_[center][i] - extrap_tab_[up][i]) * BS_.coef(center);
            }
        }

        Scalar calc_error_of_row(size_t k) const {
            size_t center = at(k, k);
            size_t left = center - 1;
            Scalar max_err = 0;

            for (size_t i = 0; i < var_num_; ++i) {
                Scalar d = space::abs(extrap_tab_[center][i] - extrap_tab_[left][i]);
                Scalar scale = space::min(space::abs(extrap_tab_[left][i]),
                                          space::abs(extrap_tab_[center][i])) + abs_error_;
                max_err = space::max(1.0 * max_err, d / scale);
            }
            return max_err / rel_error_;
        }

        Scalar calc_ideal_step_coef(Scalar error, size_t k) {
            if (error != 0) {
                return 0.9 * pow(0.9 / error, BS_.expon(k));
            } else {
                return 1.0 / BS_.limiter(k);
            }
            /*return SpaceH::max(BS_.limiter(k) / 4,
                               SpaceH::min(0.9 * pow(0.90 / error, BS_.expon(k)), 1.0 / BS_.limiter(k)));*/
        }

        inline bool in_converged_window(size_t iter) {
            return iter == ideal_iter_ - 1 || iter == ideal_iter_ || iter == ideal_iter_ + 1;
        }

        inline bool is_diverged_anyhow(Scalar error, size_t iter) const {
            Scalar r = 1;

            if (iter == ideal_iter_ - 1) {
                r = static_cast<Scalar>(BS_.step(iter + 1) * BS_.step(iter + 2)) /
                    static_cast<Scalar>(BS_.step(0) * BS_.step(0));
            } else if (iter == ideal_iter_) {
                r = static_cast<Scalar>(BS_.step(iter + 1)) / static_cast<Scalar>(BS_.step(0));
            } else {
                return error > 1.0;//k == iterDepth+1 and error >1 reject directly
            }

            return error > r * r;
        }

        inline Scalar step_coef_limiter(Scalar step_coef) {
            return space::in_range(static_cast<Scalar>(BS_.limiter(ideal_iter_) / 4), step_coef,
                                   static_cast<Scalar>(1.0 / BS_.limiter(ideal_iter_)));
        }

        inline Scalar reduced_step_coef(size_t iter) {
            if (iter == ideal_iter_ - 1)
                return optm_step_coef_[iter] * static_cast<Scalar>(BS_.work(iter + 1)) /
                       static_cast<Scalar>(BS_.work(iter));
            else
                return optm_step_coef_[ideal_iter_];
        }

        inline size_t allowed(size_t i) const {
            return space::in_range(static_cast<size_t>(2), i, static_cast<size_t>(max_depth_ - 1));
        }

        Scalar prepare_next_iter(size_t iter) {
            switch (static_cast<int>(iter - ideal_iter_)) {
                case -1:
                    //ideal_iter_ <= 2 here to avoid none calculated work_per_len[1-1=0]
                    if (work_per_len_[iter] < 0.9 * work_per_len_[iter - 1] || ideal_iter_ <= 2) {
                        return optm_step_coef_[iter] * static_cast<Scalar>(BS_.work(iter + 1)) /
                               static_cast<Scalar>(BS_.work(iter));
                    } else {
                        ideal_iter_ = allowed(ideal_iter_ - 1);
                        return optm_step_coef_[ideal_iter_];//reduce order
                    }

                case 0:
                    if (work_per_len_[iter - 1] < 0.8 * work_per_len_[iter]) {
                        ideal_iter_ = allowed(ideal_iter_ - 1);
                        return optm_step_coef_[ideal_iter_];
                    } else if (work_per_len_[iter] < 0.9 * work_per_len_[iter - 1]) {
                        ideal_iter_ = allowed(ideal_iter_ + 1);
                        return optm_step_coef_[iter] * static_cast<Scalar>(BS_.work(iter + 1)) /
                               static_cast<Scalar>(BS_.work(iter));
                    } else {
                        return optm_step_coef_[ideal_iter_];//keep order
                    }

                case 1:
                    if (work_per_len_[iter - 2] < 0.8 * work_per_len_[iter - 1]) {
                        ideal_iter_ = allowed(ideal_iter_ - 1);
                    }
                    if (work_per_len_[iter] < 0.9 * work_per_len_[ideal_iter_]) {
                        ideal_iter_ = allowed(ideal_iter_ + 1);
                    }
                    return optm_step_coef_[ideal_iter_];

                default: SPACEHUB_ABORT("unexpected iteration index!");
            }
        }

    private:
        /** @brief The Maximum iteration depth*/
        static constexpr size_t max_depth_{8};

        /** @brief The constat coef for BS extrapolation*/
        BStab<Scalar, max_depth_ + 1> BS_;

        /** @brief Extrapolation table.*/
        std::array<std::vector<Scalar>, (max_depth_ + 1) * (max_depth_ + 2) / 2> extrap_tab_;

        /** @brief The optimal step size array.*/
        std::array<Scalar, max_depth_ + 1> optm_step_coef_;

        /** @brief The work(computation resource) needed to converge at column k.*/
        std::array<Scalar, max_depth_ + 1> work_per_len_;

        /** @brief Local absolute error*/
        Scalar abs_error_{max(space::epsilon_v<Scalar>, 1e-13)};

        /** @brief Local relative error*/
        Scalar rel_error_{max(space::epsilon_v<Scalar>, 1e-13)};

        /** @brief Current iteraation depth.*/
        size_t ideal_iter_{7};

        /** @brief Total volume of extrapolation table(in scalar).*/
        size_t var_num_{0};

        /** @brief Rejectin number*/
        size_t rej_num_{0};

        /** @brief Total iteration number*/
        size_t iter_num_{0};
    };
}
#endif //SPACEHUB_BS_ITERATOR0_H
