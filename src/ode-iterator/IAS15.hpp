
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
        Scalar iterate(U &particles, Scalar macro_step_size);

        void set_rtol(Scalar rtol) { PC_err_checker_.set_rtol(rtol); };

       private:
        inline void reset_PC_iteration();

        bool in_converged_window();

        Integrator integrator_;
        StepController step_controller_;
        ErrEstimator err_checker_;
        ErrEstimator PC_err_checker_;
        Scalar last_PC_error_{math::max_value<Scalar>::value};
        Scalar last_error_{1};
        static constexpr size_t max_iter_{30};
        bool warmed_up{false};
    };

    /*---------------------------------------------------------------------------*\
          Class IAS15 Implementation
    \*---------------------------------------------------------------------------*/
    template <typename Integrator, typename ErrEstimator, typename StepController>
    IAS15<Integrator, ErrEstimator, StepController>::IAS15() {
        PC_err_checker_.set_atol(0);
        PC_err_checker_.set_rtol(1e-16);
        err_checker_.set_atol(0);
        err_checker_.set_rtol(5e-10);
        step_controller_.set_safe_guards(0.85, 1, 6e-5, 1);  //(1/6e-5)**(1/7)~4
    }

    template <typename Integrator, typename ErrEstimator, typename StepController>
    template <typename U>
    auto IAS15<Integrator, ErrEstimator, StepController>::iterate(U &particles, Scalar macro_step_size) -> Scalar {
        Scalar iter_h = macro_step_size;
        for (size_t k = 0; k < max_iter_; ++k) {
            integrator_.correct(particles, iter_h);
            if (in_converged_window()) {
                Scalar error = err_checker_.error(integrator_.y_h(), integrator_.b()[6]);

                Scalar new_iter_h = step_controller_.next_step_size((Integrator::order - 1) / 2, iter_h, error);

                if (error < 1) {
                    integrator_.evaluate(particles, iter_h);
                    integrator_.predict(new_iter_h / iter_h);
                    last_error_ = error;
                    warmed_up = true;
                    return new_iter_h;
                } else {
                    if (warmed_up) {
                        integrator_.predict(new_iter_h / iter_h);
                    }
                    iter_h = new_iter_h;
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
        Scalar PC_error = PC_err_checker_.error(integrator_.y_h(), integrator_.diff_b6());
        if (PC_error < static_cast<Scalar>(1) || PC_error >= last_PC_error_) {
            reset_PC_iteration();
            return true;
        } else {
            last_PC_error_ = PC_error;
            return false;
        }
    }
}  // namespace space::ode_iterator
