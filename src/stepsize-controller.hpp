//
// Created by 王艺涵 on 10/9/19.
//

#ifndef SPACEHUB_STEPSIZE_CONTROLLER_HPP
#define SPACEHUB_STEPSIZE_CONTROLLER_HPP


#include "dev-tools.hpp"

namespace space {

/*---------------------------------------------------------------------------*\
    Class StepController Declaration
\*---------------------------------------------------------------------------*/
/**
 * Abstract class of step size controller. A class implements(partly/fully) the interfaces of this
 * class via CRTP idiom can be used cross the system as an implementation of the concept `StepController`. The step size
 * controller provides the interface to estimate the next step size with current step size and errors in previous steps.
 *
 * @tparam Derived The implement class in CRTP idiom.
 */
  template<typename Derived>
  class StepController {
  public:
    // public methods

    /**
     * @auto_impl
     *
     * The downcast interface of Base class to Derived class.
     * @return Derived
     */
    Derived &derived();


    /**
     * @must_impl
     *
     * Estimate the step size of the next step from errors in previous steps and current step size of an integration system.
     *
     * @tparam Scalar Floating point like type type.
     * @tparam Array Tuple/array like type. i.e `std::tuple`, `std::array`.
     * @param[in] order Order of the integration method.
     * @param[in] old_step The current step size of the integration system.
     * @param[in] errors Errors of previous steps.
     * @return The estimated step size of the next step.
     */
    template<typename Scalar, typename Array>
    Scalar next_step_size(size_t order, Scalar old_step, Array const &errors);

  private:
    /**
     * Construct a new StepController object
     */
    StepController() = default;

    friend Derived;
  };

/*---------------------------------------------------------------------------*\
    Class Particles Implementation
\*---------------------------------------------------------------------------*/
  template<typename Derived>
  Derived &StepController<Derived>::derived() {
    return static_cast<Derived &>(*this);
  }

  template<typename Derived>
  template<typename Scalar, typename Array>
  Scalar StepController<Derived>::next_step_size(size_t order, Scalar old_step, Array const &errors) {
    return static_cast<Derived *>(this)->impl_next_step_size(order, old_step, errors);
  }
}
#endif //SPACEHUB_STEPSIZE_CONTROLLER_HPP
