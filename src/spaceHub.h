#ifndef SPACEHUB_H
#define SPACEHUB_H

#include "kahan-number.hpp"
#include "macros.hpp"

#include "particles/finite-size.hpp"
#include "particles/point-particles.hpp"

#include "simulator.hpp"

#include "interaction/newtonian.hpp"

#include "particle-system/archain.hpp"
#include "particle-system/base-system.hpp"
#include "particle-system/chain-system.hpp"
#include "particle-system/regu-system.hpp"

#include "integrator/symplectic/symplectic-integrator.hpp"

#include "error-checker/worst-offender.hpp"
#include "error-checker/RMS.hpp"

#include "step-controller/PID-controller.hpp"

#include "ode-iterator/Burlish-Stoer.hpp"
#include "ode-iterator/const-iterator.hpp"

#include "args-callback/callbacks.hpp"

#include "orbits/orbits.hpp"

#include "tools/auto-name.hpp"

namespace space {
  using DefaultTypes = Types<double, std::vector>;

  template<template<class> class Paticles = PointParticles, typename Force = interactions::NewtonianGrav>
  using DefaultSolver =
  Simulator<ARchainSystem<Paticles<DefaultTypes>, Force, ReguType::logH>, odeIterator::BurlishStoer<double, RMS, PIDController>>;

// template<template <class> class Paticles = SoAPointParticles, typename Force = interactions::NewtonianGrav>
// using DefaultSolver = Simulator<ChainSystem<Paticles<DefaultTypes>, Force>, odeIterator::BurlishStoer<double>>;

// template<template <class> class Paticles = SoAPointParticles, typename Force = interactions::NewtonianGrav>
// using DefaultSolver = Simulator<SimpleSystem<Paticles<DefaultTypes>, Force>, odeIterator::BurlishStoer<double>>;
}  // namespace space

#endif
