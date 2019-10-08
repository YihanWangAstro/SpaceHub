#ifndef BSITERATOR_H
#define BSITERATOR_H

#include <array>
#include <vector>
#include "ode-iterator.hpp"
#include "../own-math.hpp"

namespace space::odeIterator {

    /*---------------------------------------------------------------------------*\
        Class BSIterator Declaration
    \*---------------------------------------------------------------------------*/
    template<typename Real>
    class BSIterator : public OdeIterator<BSIterator<Real>> {
        static_assert(std::is_floating_point<Real>::value, "Only float-like type can be used!");
    public:
        //Type members
        using Scalar = Real;

        using Base = OdeIterator<BSIterator<Real>>;

        /*---------------------------------------------------------------------------*\
        Sub-Class Burlish_Stoer_tab Declaration
        \*---------------------------------------------------------------------------*/
        template<size_t MaxIter>
        class BStab {
        public:
            //Constructors
            BStab();

            //Public members
            /** @brief Extrapolation coefficient.*/
            std::array<Scalar, MaxIter * (MaxIter + 1) / 2> coef;

            /** @brief The exponent of error estimate at column k.*/
            std::array<Scalar, MaxIter> expon;

            /** @brief The minimal coeffient of integration step estimation.*/
            std::array<Scalar, MaxIter> min_step_coef;

            /** @brief The maximal coeffient of integration step estimation.*/
            std::array<Scalar, MaxIter> max_step_coef;

            /** @brief The work(computation resource) per step size of each iteration depth.*/
            std::array<size_t, MaxIter> work;

            std::array<Scalar, MaxIter> work_ratio;

            std::array<Scalar, MaxIter> step_ratio;

            /** @brief Steps of integration of each iteration depth.*/
            std::array<size_t, MaxIter> step;

        private:
            //Private static members
            static constexpr Scalar S3{0.02};

            static constexpr Scalar S4{4.0};

            //Private methods
            size_t at(size_t i, size_t j) const;
        };

        //Public methods
        void set_abs_err(Scalar err);

        void set_rel_err(Scalar err);

    CRTP_impl:
        //CRTP implementation
        template<typename U>
        Scalar impl_iterate(U &ptcs, Scalar macro_step_size);

    private:
        //Private methods
        template<typename U>
        void evolve_by_n_steps(U &ptcs, Scalar macro_step_size, size_t steps);

        void extrapolate_tab(size_t k);

        Scalar calc_error_of_row(size_t k) const;

        Scalar calc_ideal_step_coef(Scalar error, size_t k);

        Scalar step_coef_limiter(Scalar step_coef);

        Scalar ideal_step_coef(size_t ideal_k, size_t k);

        Scalar prepare_next_iter(size_t iter);

        size_t at(size_t i, size_t j) const;

        size_t allowed(size_t i) const;

        bool is_diverged_anyhow(Scalar error, size_t iter) const;

        bool in_converged_window(size_t iter);

    private:
        //Static private members
        static constexpr Scalar S1{0.95};

        static constexpr Scalar S2{0.95};

        /** @brief The Maximum iteration depth*/
        static constexpr size_t max_depth{8};

        static constexpr bool most_offensive_err{true};

        //Private members
        /** @brief The constat coef for BS extrapolation*/
        BStab<max_depth + 1> BS_;

        /** @brief Extrapolation table.*/
        std::array<std::vector<Scalar>, (max_depth + 1) * (max_depth + 2) / 2> extrap_tab_;

        /** @brief The optimal step size array.*/
        std::array<Scalar, max_depth + 1> optm_step_coef_;

        /** @brief The work(computation resource) needed to converge at column k.*/
        std::array<Scalar, max_depth + 1> work_per_len_;

        /** @brief Local absolute error*/
        Scalar abs_error_{max(space::epsilon_v<Scalar>, 1e-13)};

        /** @brief Local relative error*/
        Scalar rel_error_{max(space::epsilon_v<Scalar>, 1e-14)};

        /** @brief Current iteraation depth.*/
        size_t ideal_iter_{7};

        /** @brief Total volume of extrapolation table(in scalar).*/
        size_t var_num_{0};

        /** @brief Rejectin number*/
        size_t rej_num_{0};

        /** @brief Total iteration number*/
        size_t iter_num_{0};
    };

    /*---------------------------------------------------------------------------*\
        Class BSIterator Implementation
    \*---------------------------------------------------------------------------*/
    template<typename Real>
    void BSIterator<Real>::set_abs_err(Scalar err) {
        abs_error_ = err;
    }

    template<typename Real>
    void BSIterator<Real>::set_rel_err(Scalar err) {
        rel_error_ = err;
    }

    template<typename Real>
    template<typename U>
    auto BSIterator<Real>::impl_iterate(U &ptcs, Scalar macro_step_size) -> Scalar {
        static_assert(is_particle_system_v<U>, "Passing non paritcle-system-type!");

        Scalar iter_H = macro_step_size;

        for (;;) {
            iter_num_++;
            auto local_sys = ptcs;

            evolve_by_n_steps(local_sys, iter_H, BS_.step[0]);

            local_sys.to_linear_container(extrap_tab_[at(0, 0)]);

            var_num_ = extrap_tab_[0].size();

            for (size_t i = 1; i <= ideal_iter_ + 1; ++i) {
                local_sys = ptcs;
                evolve_by_n_steps(local_sys, iter_H, BS_.step[i]);
                local_sys.to_linear_container(extrap_tab_[at(i, 0)]);

                extrapolate_tab(i);

                Scalar error = calc_error_of_row(i);
                optm_step_coef_[i] = calc_ideal_step_coef(error, i);
                work_per_len_[i] = BS_.work[i] / optm_step_coef_[i];

                if (in_converged_window(i)) {
                    if (error < 1.0) {
                        iter_H *= step_coef_limiter(prepare_next_iter(i));
                        ptcs.load_from_linear_container(extrap_tab_[at(i, i)]);
                        return iter_H;
                    } else if (is_diverged_anyhow(error, i)) {
                        rej_num_++;
                        iter_H *= step_coef_limiter(ideal_step_coef(ideal_iter_, i));
                        break;
                    }
                }
            }
        }
    }

    template<typename Real>
    template<typename U>
    void BSIterator<Real>::evolve_by_n_steps(U &ptcs, Scalar macro_step_size, size_t steps) {
        Scalar h = macro_step_size / steps;

        ptcs.drift(0.5 * h);
        for (size_t i = 1; i < steps; i++) {
            ptcs.kick(h);
            ptcs.drift(h);
        }
        ptcs.kick(h);
        ptcs.drift(0.5 * h);
    }

    template<typename Real>
    void BSIterator<Real>::extrapolate_tab(size_t k) {
        size_t curr_row = at(k, 0);
        size_t last_row = at(k - 1, 0);

        for (size_t j = 0; j < k; ++j) {
            size_t center = curr_row + j;//
            size_t right = center + 1;
            size_t up = last_row + j;

            if (extrap_tab_[right].size() != var_num_)
                extrap_tab_[right].resize(var_num_);

            for (size_t i = 0; i < var_num_; ++i)
                extrap_tab_[right][i] =
                        extrap_tab_[center][i] + (extrap_tab_[center][i] - extrap_tab_[up][i]) * BS_.coef[center];
        }
    }

    template<typename Real>
    auto BSIterator<Real>::calc_error_of_row(size_t k)  const -> Scalar {
        size_t center = at(k, k);
        size_t left = center - 1;
        Scalar max_err = 0;

        if constexpr (most_offensive_err){
            for (size_t i = 0; i < var_num_; ++i) {
                Scalar d = space::abs(extrap_tab_[center][i] - extrap_tab_[left][i]);
                Scalar y = space::min(space::abs(extrap_tab_[left][i]), space::abs(extrap_tab_[center][i]));
                if(iseq(y, static_cast<Scalar>(0)))
                    continue;
                Scalar scale = rel_error_* space::min(space::abs(extrap_tab_[left][i]), space::abs(extrap_tab_[center][i]));
                max_err = space::max(1.0 * max_err, d / scale);
            }
        } else {
            for (size_t i = 0; i < var_num_; ++i) {
                Scalar d = extrap_tab_[center][i] - extrap_tab_[left][i];
                Scalar scale = rel_error_* space::min(space::abs(extrap_tab_[left][i]), space::abs(extrap_tab_[center][i])) + abs_error_;
                max_err += d*d/(scale*scale);
            }

            max_err = sqrt(max_err/var_num_);
        }

        return max_err;
    }

    template<typename Real>
    auto BSIterator<Real>::calc_ideal_step_coef(Real error, size_t k) -> Scalar {
        if (error != 0) {
            return S1 * pow(S2 / error, BS_.expon[k]);
        } else {
            return BS_.max_step_coef[k];
        }
        /*return SpaceH::max(BS_.limiter(k) / 4,
                           SpaceH::min(0.9 * pow(0.90 / error, BS_.expon(k)), 1.0 / BS_.step_limiter(k)));*/
    }

    template<typename Real>
    bool BSIterator<Real>::is_diverged_anyhow(Real error, size_t iter) const {
        Real r = 1;

        if (iter == ideal_iter_ - 1) {
            r = BS_.step_ratio[iter + 1] * BS_.step_ratio[iter + 2];
        } else if (iter == ideal_iter_) {
            r = BS_.step_ratio[iter + 1];
        } else {
            return error > 1.0;//k == iterDepth+1 and error >1 reject directly
        }

        return error > r * r;
    }

    template<typename Real>
    auto BSIterator<Real>::ideal_step_coef(size_t ideal_k, size_t k) -> Scalar {
        if(ideal_k == k + 1){
            return optm_step_coef_[k] * BS_.work_ratio[k + 1];
        } else if (ideal_k <= k) {
            return optm_step_coef_[ideal_k];
        } else {
            spacehub_abort("Unexpected iteration steps!");
            return 0;
        }
    }

    template<typename Real>
    auto BSIterator<Real>::prepare_next_iter(size_t iter) -> Scalar {
        if(iter == ideal_iter_ - 1 || iter == ideal_iter_ ) {
            if(work_per_len_[iter - 1] < 0.8* work_per_len_[iter]){
                ideal_iter_ = allowed(iter - 1);
            } else if(work_per_len_[iter] < 0.9 * work_per_len_[iter - 1]) {
                ideal_iter_ = allowed(iter + 1);
            } else {
                ideal_iter_ = allowed(iter);
            }
        } else {
            if(work_per_len_[iter - 2] < 0.8 * work_per_len_[iter - 1]) {
                ideal_iter_ = allowed(ideal_iter_ - 1);
            }
            if (work_per_len_[iter] < 0.9* work_per_len_[ideal_iter_]) {
                ideal_iter_ = allowed(ideal_iter_ + 1);
            }
        }
        return ideal_step_coef(ideal_iter_, iter);
    }

    template<typename Real>
    inline bool BSIterator<Real>::in_converged_window(size_t iter) {
        return iter == ideal_iter_ - 1 || iter == ideal_iter_ || iter == ideal_iter_ + 1;
    }

    template<typename Real>
    inline auto BSIterator<Real>::step_coef_limiter(Scalar step_coef) -> Scalar {
        return space::in_range(BS_.min_step_coef[ideal_iter_], step_coef, BS_.max_step_coef[ideal_iter_]);
    }

    template<typename Real>
    inline size_t BSIterator<Real>::allowed(size_t i) const {
        return space::in_range(static_cast<size_t>(2), i, static_cast<size_t>(max_depth - 1));
    }

    template<typename Real>
    inline size_t BSIterator<Real>::at(size_t i, size_t j) const {
        return i * (i + 1) / 2 + j;
    }

    /*---------------------------------------------------------------------------*\
       Sub-Class Burlish_Stoer_tab Implementation
    \*---------------------------------------------------------------------------*/
    template<typename Real>
    template<size_t MaxIter>
    BSIterator<Real>::BStab<MaxIter>::BStab() {
        for (size_t i = 0; i < MaxIter; ++i) {

            step[i] = 2 * (i + 1);

            step_ratio[i] = static_cast<Scalar>(step[i])/ static_cast<Scalar>(step[0]);

            work[i] = (i==0 ?  0 : work[i - 1]) + step[i];

            work_ratio[i] = static_cast<Scalar>(work[i]) / static_cast<Scalar>(i == 0 ? work[i] : work[i-1]);

            expon[i] = static_cast<Scalar>(1.0 / (2 * i + 1));

            auto F = pow(S3, expon[i]);

            min_step_coef[i] = F / S4;
            max_step_coef[i] = 1.0 / F;

            for (size_t j = 0; j < i; ++j) {
                Scalar ratio = static_cast<Scalar>(step[i]) / static_cast<Scalar>(step[i - j - 1]);
                coef[at(i, j)] = 1.0 / (ratio * ratio - 1);
            }
        }
    }

    template<typename Real>
    template<size_t MaxIter>
    inline size_t BSIterator<Real>::BStab<MaxIter>::at(size_t i, size_t j) const {
        return i * (i + 1) / 2 + j;
    }
}

#endif
