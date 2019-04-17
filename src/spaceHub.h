#ifndef SPACEHUB_H
#define SPACEHUB_H

#include "macros.h"
#include "kahan-number.h"

#include "particle/point-particle.tpp"
#include "particle/finite-size.tpp"

#include "simulator.h"

#include "interaction/newtonian.h"

#include "particle-system/base-system.tpp"
#include "particle-system/regu-system.tpp"
#include "particle-system/chain-system.tpp"
#include "particle-system/archain.tpp"

#include "integrator/symplectic/symplectic-integrator.tpp"

#include "ode-iterator/const-iterator.tpp"
#include "ode-iterator/BS-iterator.tpp"

#include "args-callback/callbacks.h"

#include "orbits/orbits.h"

#include "tools/auto-name.tpp"

namespace space{
    using DefaultTypes = Types<double, std::vector>;

    template<template <class> class Paticles, typename Force>
    using DefaultSolver = Simulator<ARchainSystem<Paticles<DefaultTypes>, Force, ReguType::logH>, odeIterator::BSIterator<double>>;
}

#endif
