
#ifndef SPACEHUB_CONST_ITERATOR_HPP
#define SPACEHUB_CONST_ITERATOR_HPP

#include "ode-iterator.hpp"

namespace space::odeIterator {

  /*---------------------------------------------------------------------------*\
      Class ConstOdeIterator Declaration
  \*---------------------------------------------------------------------------*/
  template<typename Integrator>
  class ConstOdeIterator : public OdeIterator<ConstOdeIterator<Integrator>> {
  public:
    //Type members
    using Base = OdeIterator<ConstOdeIterator<Integrator>>;

    CRTP_IMPL:

    //CRTP implementation
    template<typename T>
    auto impl_iterate(T &particles, typename T::Scalar macro_step_size) -> typename T::Scalar;

  private:
    //Private members
    Integrator integrator_;
  };

  /*---------------------------------------------------------------------------*\
      Class ConstOdeIterator Implementation
  \*---------------------------------------------------------------------------*/
  template<typename Integrator>
  template<typename T>
  auto
  ConstOdeIterator<Integrator>::impl_iterate(T &particles, typename T::Scalar macro_step_size) -> typename T::Scalar {
    static_assert(is_particle_system_v<T>, "Passing non particle-system-type!");
    integrator_.integrate(particles, macro_step_size);
    return macro_step_size;
  }
}
#endif
