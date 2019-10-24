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
    the terms of the MIT License. SpaceHub is distributed in the hope that it
    will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the MIT License
    for more details. You should have received a copy of the MIT License along
    with SpaceHub.
\*---------------------------------------------------------------------------*/
/**
 * @file const-iterator.hpp
 *
 * Header file.
 */
#ifndef SPACEHUB_CONST_ITERATOR_HPP
#define SPACEHUB_CONST_ITERATOR_HPP

#include "ode-iterator.hpp"

namespace space::ode_iterator {

/*---------------------------------------------------------------------------*\
      Class ConstOdeIterator Declaration
\*---------------------------------------------------------------------------*/
  /**
   *
   * @tparam Integrator
   */
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
    static_assert(particle_system::is_particle_system_v<T>, "Passing non particle-system-type!");
    integrator_.integrate(particles, macro_step_size);
    return macro_step_size;
  }
}
#endif
