#ifndef SPACEHUB_H
#define SPACEHUB_H

#include "macros.hpp"
#include "kahan-number.hpp"

#include "particle/point-particles.hpp"
#include "particle/finite-size.hpp"

#include "simulator.hpp"

#include "interaction/newtonian.hpp"

#include "particle-system/base-system.hpp"
#include "particle-system/regu-system.hpp"
#include "particle-system/chain-system.hpp"
#include "particle-system/archain.hpp"

#include "integrator/symplectic/symplectic-integrator.hpp"

#include "ode-iterator/const-iterator.hpp"
#include "ode-iterator/BS-iterator.hpp"

#include "args-callback/callbacks.hpp"

#include "orbits/orbits.hpp"

#include "tools/auto-name.hpp"

namespace space{
    using DefaultTypes = Types<double, std::vector>;

    template<template <class> class Paticles = PointParticles, typename Force = interactions::NewtonianGrav>
    using DefaultSolver = Simulator<ARchainSystem<Paticles<DefaultTypes>, Force, ReguType::TTL>, odeIterator::BSIterator<double>>;

    //template<template <class> class Paticles = SoAPointParticles, typename Force = interactions::NewtonianGrav>
    //using DefaultSolver = Simulator<ChainSystem<Paticles<DefaultTypes>, Force>, odeIterator::BSIterator<double>>;

    //template<template <class> class Paticles = SoAPointParticles, typename Force = interactions::NewtonianGrav>
    //using DefaultSolver = Simulator<SimpleSystem<Paticles<DefaultTypes>, Force>, odeIterator::BSIterator<double>>;
}

#endif
