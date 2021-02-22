
#pragma once
/**
 * @file IAS15.hpp
 *
 * Header file.
 */
#include "../dev-tools.hpp"
#include "../integrator/Gauss-Radau.hpp"
#include "../math.hpp"

namespace space::ode_iterator {
    /*---------------------------------------------------------------------------*\
          Class IAS15 Declaration
    \*---------------------------------------------------------------------------*/

    /**
     * @brief IAS15 iterator see details in https://arxiv.org/abs/1409.4779 .
     *
     * @tparam Integrator
     * @tparam ErrEstimator
     * @tparam StepController
     */
    template <typename Integrator, typename ErrEstimator, typename StepController>
    class IAS15 {
       public:
        SPACEHUB_USING_TYPE_SYSTEM_OF(Integrator);
        static_assert(std::is_same_v<Integrator, integrator::GaussRadau<TypeSet>>,
                      "IAS15 iterator only works with Gauss-Radau integrator!");

        IAS15();

        template <typename U>
        Scalar iterate(U& particles, Scalar macro_step_size);

        void set_rtol(Scalar rtol) { PC_err_checker_.set_rtol(rtol); };

       private:
        inline void reset_PC_iteration();

        bool in_converged_window();

        template <typename Array1, typename Array2, typename U>
        Scalar calc_step_error(Array1 const& dy_h, Array2 const& b6, U const& ptc, Scalar step_size) const;

        Integrator integrator_;
        StepController step_ctrl_;
        ErrEstimator PC_err_checker_;
        Scalar last_PC_error_{math::max_value<Scalar>::value};
        static constexpr size_t max_iter_{30};
        bool warmed_up{false};

        CREATE_STATIC_MEMBER_CHECK(regu_type);
        CREATE_METHOD_CHECK(chain_pos);
        CREATE_METHOD_CHECK(chain_vel);
    };

    /*---------------------------------------------------------------------------*\
          Class IAS15 Implementation
    \*---------------------------------------------------------------------------*/
    template <typename Integrator, typename ErrEstimator, typename StepController>
    IAS15<Integrator, ErrEstimator, StepController>::IAS15() {
        PC_err_checker_.set_atol(0);
        PC_err_checker_.set_rtol(1e-16);
        step_ctrl_.set_safe_guards(0.85, 1.0);
        step_ctrl_.set_limiter(0.02, 4.0);
    }

    template <typename Integrator, typename ErrEstimator, typename StepController>
    template <typename U>
    auto IAS15<Integrator, ErrEstimator, StepController>::iterate(U& particles, Scalar macro_step_size) -> Scalar {
        Scalar iter_h = macro_step_size;
        for (size_t k = 0; k < max_iter_; ++k) {
            integrator_.correct(particles, iter_h);
            if (in_converged_window()) {
                Scalar step_error = calc_step_error(integrator_.dy_h(), integrator_.b()[6], particles, iter_h);
                Scalar new_step_ratio = step_ctrl_.next((Integrator::order - 1) / 2, step_error);
                // print_csv(std::cout, k, error, iter_h, '\n');
                if (new_step_ratio > step_ctrl_.limiter_min()) {
                    integrator_.evaluate(particles, iter_h);
                    new_step_ratio = step_ctrl_.limiter(new_step_ratio);
                    integrator_.predict(new_step_ratio);
                    warmed_up = true;
                    return iter_h * new_step_ratio;
                } else {
                    if (warmed_up) {
                        integrator_.predict(new_step_ratio);
                    }
                    iter_h *= new_step_ratio;
                    k = 0;
                    reset_PC_iteration();

                    continue;
                }
            }
        }
        spacehub_abort("Exceed the max iteration number");
    }

    template <typename Integrator, typename ErrEstimator, typename StepController>
    void IAS15<Integrator, ErrEstimator, StepController>::reset_PC_iteration() {
        last_PC_error_ = math::max_value<Scalar>::value;
    }

    template <typename Integrator, typename ErrEstimator, typename StepController>
    bool IAS15<Integrator, ErrEstimator, StepController>::in_converged_window() {
        Scalar PC_error = PC_err_checker_.error(integrator_.dy_h(), integrator_.diff_b6());
        if (PC_error < static_cast<Scalar>(1) || PC_error >= last_PC_error_) {
            reset_PC_iteration();
            return true;
        } else {
            last_PC_error_ = PC_error;
            return false;
        }
    }

    template <typename Integrator, typename ErrEstimator, typename StepController>
    template <typename Array1, typename Array2, typename U>
    auto IAS15<Integrator, ErrEstimator, StepController>::calc_step_error(Array1 const& dy_h, Array2 const& b6,
                                                                          U const& ptc, Scalar step_size) const
        -> Scalar {
        static Scalar step_rtol_{5e-10};
        Scalar max_diff = 0;
        Scalar max_scale = 0;
        size_t size = dy_h.size();
        size_t ptc_num = ptc.number();
        size_t pos_offset = ptc.pos_offset();
        size_t vel_offset = ptc.vel_offset();
        size_t auxi_vel_offset = ptc.auxi_vel_offset();

        Scalar dt = step_size / ptc.step_scale();
        Scalar dt2 = dt * dt;
        std::vector<bool> mask(size, false);
        for (size_t i = 0; i < ptc_num; i++) {
            bool slow_varing = false;
            if constexpr (HAS_METHOD(U, chain_pos) && HAS_METHOD(U, chain_vel)) {
                slow_varing = norm2(ptc.chain_pos(i)) * 1e-12 > norm2(ptc.chain_vel(i)) * dt2;
            } else {
                slow_varing = norm2(ptc.pos(i)) * 1e-12 > norm2(ptc.vel(i)) * dt2;
            }
            if (norm2(ptc.pos(i)) * 1e-12 > norm2(ptc.vel(i)) * dt2) {
                mask[pos_offset + 3 * i] = mask[pos_offset + 3 * i + 1] = mask[pos_offset + 3 * i + 2] = true;
                mask[vel_offset + 3 * i] = mask[vel_offset + 3 * i + 1] = mask[vel_offset + 3 * i + 2] = true;
                if constexpr (U::ext_vel_dep) {
                    mask[auxi_vel_offset + 3 * i] = mask[auxi_vel_offset + 3 * i + 1] =
                        mask[auxi_vel_offset + 3 * i + 2] = true;
                }
            }
        }

        for (size_t i = 0; i < size; i++) {
            if (mask[i]) {
                continue;
            }
            max_diff = std::max(max_diff, static_cast<Scalar>(fabs(b6[i])));
            max_scale = std::max(max_scale, static_cast<Scalar>(fabs(dy_h[i])));
        }

        if (max_scale == 0) {
            return 0;
        } else {
            return max_diff / (max_scale * step_rtol_);
        }
    }
}  // namespace space::ode_iterator
