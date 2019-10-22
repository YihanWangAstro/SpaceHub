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
 * @file integrator.hpp
 *
 * Header file.
 */
#ifndef SPACEHUB_INTEGRATOR_HPP
#define SPACEHUB_INTEGRATOR_HPP

#include "../particle-system/particle-system.hpp"

namespace space::integrator {
  /*---------------------------------------------------------------------------*\
    Class Integrator Declaration
\*---------------------------------------------------------------------------*/

  /**
   * Abstract class of integrator. A class implements(partly/fully) the interfaces of this
   * class via CRTP idiom can be used cross the system as an implementation of the concept `Integrator`.
   *
   * @tparam Derived The implement class in CRTP idiom.
   */
  template<typename Derived>
  class Integrator {
  public:
    /**
     * Order of the integrator.
     */
    static constexpr size_t order{Derived::order};

    /**
     * @must_impl
     *
     * Integrate a particle system with given step_size.
     *
     * @tparam T Any implementation of concept particle_system::ParticleSystem.
     * @param[in,out] particle_sys The particle system that need integration.
     * @param[in] step_size Step size.
     */
    template<typename T>
    void integrate(T &particle_sys, typename T::Scalar step_size);

  private:
    Integrator() = default;

    friend Derived;
  };

  /*---------------------------------------------------------------------------*\
    Class Integrator Definition
\*---------------------------------------------------------------------------*/
  template<typename Derived>
  template<typename T>
  void Integrator<Derived>::integrate(T &particle_system, Scalar step_size) {
    static_assert(particle_system::is_particle_system_v<T>, "Passing non paritcle-system-type!");
    static_cast<Derived *>(this)->impl_integrate(particle_system, step_size);
  }
}
#endif //SPACEHUB_INTEGRATOR_HPP
