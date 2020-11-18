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
#pragma once
#include "../spacehub-concepts.hpp"
namespace space::ode_iterator {

/*---------------------------------------------------------------------------*\
      Class ConstOdeIterator Declaration
\*---------------------------------------------------------------------------*/
/**
 *
 * @tparam Integrator
 */
template <typename Integrator>
class ConstOdeIterator {
 public:
  template <concepts::ParticleSystem T>
  auto iterate(T &particles, typename T::Scalar macro_step_size) -> typename T::Scalar;

 private:
  // Private members
  Integrator integrator_;
};

/*---------------------------------------------------------------------------*\
      Class ConstOdeIterator Implementation
\*---------------------------------------------------------------------------*/
template <typename Integrator>
template <concepts::ParticleSystem T>
auto ConstOdeIterator<Integrator>::iterate(T &particles, typename T::Scalar macro_step_size) -> typename T::Scalar {
  integrator_.integrate(particles, macro_step_size);
  return macro_step_size;
}
}  // namespace space::ode_iterator
