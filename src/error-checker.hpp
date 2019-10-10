//
// Created by 王艺涵 on 10/8/19.
//

#ifndef SPACEHUB_ERROR_CHECKER_HPP
#define SPACEHUB_ERROR_CHECKER_HPP

#include "dev-tools.hpp"

namespace space {

/*---------------------------------------------------------------------------*\
    Class Error_checker Declaration
\*---------------------------------------------------------------------------*/
/**
 * @brief CRTP base class of a 'Structure of Array' kind particle set.
 *
 * @tparam Derived
 */
  template<typename Derived>
  class ErrorChecker {
  public:
    // public methods
    /**
     *
     * @return
     */
    DECLARE_CRTP_ACCESSOR(Derived, auto, atol);

    /**
     *
     * @return
     */
    DECLARE_CRTP_ACCESSOR(Derived, auto, rtol);

    /**
     *
     * @return
     */
    Derived &derived();

    /**
     *
     */
    template<typename T>
    void set_atol(T);

    /**
     *
     */
    template<typename T>
    void set_rtol(T);

    template<typename Array>
    auto error(Array const &y0, Array const &y1) -> typename Array::value_type;

    template<typename Array>
    auto error(Array const &scale, Array const &y0, Array const &y1) -> typename Array::value_type;

  private:
    /**
     * @brief Construct a new ErrorChecker object
     *
     */
    ErrorChecker() = default;

    friend Derived;
  };

/*---------------------------------------------------------------------------*\
    Class Particles Implementation
\*---------------------------------------------------------------------------*/
  template<typename Derived>
  Derived &ErrorChecker<Derived>::derived() {
    return static_cast<Derived &>(*this);
  }

  template<typename Derived>
  template<typename T>
  void ErrorChecker<Derived>::set_atol(T atol) {
    static_cast<Derived *>(this)->impl_set_atol(atol);
  }

  template<typename Derived>
  template<typename T>
  void ErrorChecker<Derived>::set_rtol(T rtol) {
    static_cast<Derived *>(this)->impl_set_rtol(rtol);
  }

  template<typename Derived>
  template<typename Array>
  auto ErrorChecker<Derived>::error(Array const &y0, Array const &y1) -> typename Array::value_type {
    return static_cast<Derived *>(this)->impl_error(y0, y1);
  }

  template<typename Derived>
  template<typename Array>
  auto
  ErrorChecker<Derived>::error(Array const &scale, Array const &y0, Array const &y1) -> typename Array::value_type {
    return static_cast<Derived *>(this)->impl_error(scale, y0, y1);
  }

}
#endif //SPACEHUB_ERROR_CHECKER_HPP
