//
// Created by 王艺涵 on 10/9/19.
//

#ifndef SPACEHUB_STEPSIZE_CONTROLLER_H
#define SPACEHUB_STEPSIZE_CONTROLLER_H


#include "dev-tools.hpp"

namespace space {

/*---------------------------------------------------------------------------*\
    Class StepController Declaration
\*---------------------------------------------------------------------------*/
/**
 * @brief
 *
 * @tparam Derived
 */
  template<typename Derived>
  class StepController {
  public:
    // public methods

    /**
     *
     * @return
     */
    Derived &derived();


    template<typename Scalar, typename Array>
    Scalar next_step_size(size_t order, Scalar old_step, Array const& errors);

  private:
    /**
     * @brief Construct a new StepController object
     *
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
  Scalar StepController<Derived>::next_step_size(size_t order, Scalar old_step, Array const& errors) {
    return static_cast<Derived *>(this)->impl_next_step_size(order, old_step, errors);
  }
}
#endif //SPACEHUB_STEPSIZE_CONTROLLER_H
