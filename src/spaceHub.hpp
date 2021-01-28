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
 * @dir src @brief Root
 * @folder{vector}
 * @folder{tools}
 * @folder{stellar}
 * @folder{particles}
 * @folder{particle-system}
 * @folder{orbits}
 * @folder{ode-iterator}
 *   @folder{ode-iterator/error-checker}
 *   @folder{ode-iterator/step-controller}
 * @folder{multi-thread}
 * @folder{lazy-evaluation}
 * @folder{interaction}
 * @folder{integrator}
 *   @folder{integrator/symplectic}
 * @folder{args-callback}
 * @folder{scattering}
 *
 * @file spaceHub.hpp
 *
 * Giant header file.
 */

#pragma once

#include "args-callback/callbacks.hpp"
#include "integrator/Gauss-Dadau.hpp"
#include "integrator/symplectic/symplectic-integrator.hpp"
#include "interaction/newtonian.hpp"
#include "interaction/post-newtonian.hpp"
#include "kahan-number.hpp"
#include "macros.hpp"
#include "multi-thread/multi-thread.hpp"
#include "ode-iterator/Burlish-Stoer.hpp"
#include "ode-iterator/IAS15.hpp"
#include "ode-iterator/bisec-iterator.hpp"
#include "ode-iterator/const-iterator.hpp"
#include "ode-iterator/error-checker/IAS15-error.hpp"
#include "ode-iterator/error-checker/RMS.hpp"
#include "ode-iterator/error-checker/worst-offender.hpp"
#include "ode-iterator/step-controller/PID-controller.hpp"
#include "ode-iterator/step-controller/const-controller.hpp"
#include "orbits/orbits.hpp"
#include "orbits/particle-manip.hpp"
#include "particle-system/archain.hpp"
#include "particle-system/base-system.hpp"
#include "particle-system/chain-system.hpp"
#include "particle-system/regu-system.hpp"
#include "particles/finite-size.hpp"
#include "particles/point-particles.hpp"
#include "scattering/cross-section.hpp"
#include "scattering/hierarchical.hpp"
#include "simulator.hpp"
#include "stellar/stellar.hpp"
#include "tools/auto-name.hpp"
#include "tools/config-reader.hpp"
#include "tools/timer.hpp"
#include "type-class.hpp"
/**
 * @namespace space
 * Documentation for space
 */
namespace space {

#define USING_NAMESPACE_SPACEHUB_ALL       \
    using namespace space;                 \
    using namespace space::calc;           \
    using namespace space::tools;          \
    using namespace space::ode_iterator;   \
    using namespace space::integrator;     \
    using namespace space::orbit;          \
    using namespace space::unit;           \
    using namespace space::particle_set;   \
    using namespace space::multi_thread;   \
    using namespace space::random;         \
    using namespace space::run_operations; \
    using namespace space::particle_system

    using DefaultTypes = Types<double, std::vector>;

    using DefaultSolver = Simulator<
        particle_system::ARchainSystem<particle_set::PointParticles<DefaultTypes>,
                                       interactions::Interactions<interactions::NewtonianGrav>,
                                       particle_system::ReguType::LogH>,
        ode_iterator::BurlishStoer<integrator::LeapFrogDKD<DefaultTypes>, ode_iterator::WorstOffender<DefaultTypes>,
                                   ode_iterator::PIDController<DefaultTypes>>>;
}  // namespace space
