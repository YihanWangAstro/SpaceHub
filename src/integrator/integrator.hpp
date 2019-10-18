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
  /**
   *
   * @tparam Derived
   */
  template<typename Derived>
  class Integrator {
  public:
    static constexpr size_t order{Derived::order};

    template<typename T>
    void integrate(T &particle_system, typename T::Scalar step_size) {
      static_assert(particle_system::is_particle_system_v<T>, "Passing non paritcle-system-type!");
      static_cast<Derived *>(this)->impl_integrate(particle_system, step_size);
    }

  private:
    Integrator() = default;

    friend Derived;
  };
}
#endif //SPACEHUB_INTEGRATOR_HPP
