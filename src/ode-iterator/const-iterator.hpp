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
    the terms of the GPL-3.0 License. SpaceHub is distributed in the hope that it
    will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GPL-3.0 License
    for more details. You should have received a copy of the GPL-3.0 License along
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
    template <typename Integrator, typename ErrEstimator = void, typename StepController = void>
    class ConstOdeIterator {
       public:
        template <CONCEPT_PARTICLE_SYSTEM T>
        auto iterate(T &particles, typename T::Scalar macro_step_size) -> typename T::Scalar;

       private:
        // Private members
        Integrator integrator_;
    };

    /*---------------------------------------------------------------------------*\
          Class ConstOdeIterator Implementation
    \*---------------------------------------------------------------------------*/
    template <typename Integrator, typename ErrEstimator, typename StepController>
    template <CONCEPT_PARTICLE_SYSTEM T>
    auto ConstOdeIterator<Integrator, ErrEstimator, StepController>::iterate(T &particles,
                                                                             typename T::Scalar macro_step_size) ->
        typename T::Scalar {
        integrator_.integrate(particles, macro_step_size);
        return macro_step_size;
    }
}  // namespace space::ode_iterator
