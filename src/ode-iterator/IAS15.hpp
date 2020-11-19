
#pragma once
/**
 * @file IAS15.hpp
 *
 * Header file.
 */
#include "../dev-tools.hpp"
#include "../integrator/Gauss-Dadau.hpp"
#include "../math.hpp"
#include "ode-iterator.hpp"

namespace space::ode_iterator {
    /*---------------------------------------------------------------------------*\
          Class IAS15 Declaration
    \*---------------------------------------------------------------------------*/
    /**
     * IAS15 iterator see details in https://arxiv.org/abs/1409.4779 .
     * @tparam Coords
     * @tparam ErrChecker
     * @tparam StepControl
     */
    template <typename Coords, template <typename> typename ErrChecker,
              template <size_t, typename> typename StepControl>
    class IAS15 : public OdeIterator<IAS15<Coords, ErrChecker, StepControl>> {
       public:
        using Scalar = typename Coords::Scalar;
        using Coord = Coords;
        using Integrator = integrator::GaussDadau<Coord>;

        IAS15();

        template <typename U>
        Scalar impl_iterate(U &particles, typename U::Scalar macro_step_size);

       private:
        inline void reset_PC_iteration();

        bool in_converged_window(size_t k);

        Integrator integrator_;
        typename Integrator::IterTable last_b_table_;
        StepControl<7, Scalar> step_controller_;
        ErrChecker<Scalar> err_checker_;
        ErrChecker<Scalar> PC_err_checker_;
        Scalar last_PC_error_{math::max_value<Scalar>::value};
        Scalar last_error_{1};
        static constexpr size_t max_iter_{12};
        bool warmed_up{false};
    };
    /*---------------------------------------------------------------------------*\
          Class IAS15 Implementation
    \*---------------------------------------------------------------------------*/
    template <typename Coords, template <typename> typename ErrChecker,
              template <size_t, typename> typename StepControl>
    IAS15<Coords, ErrChecker, StepControl>::IAS15() {
        PC_err_checker_.set_atol(0);
        PC_err_checker_.set_rtol(1e-16);
        err_checker_.set_atol(0);
        err_checker_.set_rtol(1e-9);
        step_controller_.set_safe_guards(0.95, 0.65, 0.02, 4);
    }

    template <typename Coords, template <typename> typename ErrChecker,
              template <size_t, typename> typename StepControl>
    template <typename U>
    auto IAS15<Coords, ErrChecker, StepControl>::impl_iterate(U &particles, typename U::Scalar macro_step_size)
        -> Scalar {
        Scalar iter_h = macro_step_size;
        integrator_.check_particle_size(particles.number());
        last_b_table_ = integrator_.b_tab();
        for (size_t k = 0; k < max_iter_; ++k) {
            integrator_.calc_B_table(particles, iter_h);
            // Scalar error = err_checker_.error(integrator_.init_acc(), integrator_.b_tab()[6]);
            // space::std_print(k, ',', iter_h, '\n');
            if (in_converged_window(k)) {
                Scalar error = err_checker_.error(integrator_.last_acc(), integrator_.b_tab()[6]);
                Scalar new_iter_h = step_controller_.next_step_size((Integrator::order - 1) / 2, iter_h,
                                                                    std::make_tuple(error, last_error_));
                // space::std_print("stp error ", k, ' ', iter_h, ' ', error, ',', new_iter_h, '\n');
                if (error < 1) {
                    integrator_.integrate_to(particles, iter_h, Integrator::final_point);
                    integrator_.predict_new_B(new_iter_h / iter_h);
                    last_error_ = error;
                    warmed_up = true;
                    return new_iter_h;
                } else {
                    if (warmed_up) {
                        integrator_.predict_new_B(new_iter_h / iter_h);
                    }
                    iter_h = new_iter_h;
                    k = 0;
                    reset_PC_iteration();
                    last_b_table_ = integrator_.b_tab();
                    continue;
                }
            }
            last_b_table_ = integrator_.b_tab();
        }
        spacehub_abort("Exceed the max iteration number");
    }

    template <typename Coords, template <typename> typename ErrChecker,
              template <size_t, typename> typename StepControl>
    void IAS15<Coords, ErrChecker, StepControl>::reset_PC_iteration() {
        last_PC_error_ = math::max_value<Scalar>::value;
    }

    template <typename Coords, template <typename> typename ErrChecker,
              template <size_t, typename> typename StepControl>
    bool IAS15<Coords, ErrChecker, StepControl>::in_converged_window(size_t k) {
        Scalar PC_error = PC_err_checker_.error(integrator_.last_acc(), last_b_table_[6], integrator_.b_tab()[6]);
        if (PC_error < static_cast<Scalar>(1) || PC_error >= last_PC_error_) {
            reset_PC_iteration();
            return true;
        } else {
            last_PC_error_ = PC_error;
            return false;
        }
    }
}  // namespace space::ode_iterator
