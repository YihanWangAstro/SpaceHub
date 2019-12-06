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
 *
 * @file collision.hpp
 *
 * Header file.
 */
#ifndef SPACEHUB_COLLISION_HPP
#define SPACEHUB_COLLISION_HPP
#include "../dev-tools.hpp"
#include "../vector/vector3.hpp"
namespace run_operations {

CREATE_METHOD_CHECK(radius);

/**
 * @brief Detect collisions with sticky radius approximation.
 *
 */
class StickyCollision {
  /**
   * @brief Callable interface for collision detection.
   *
   * @tparam ParticleSys Any implementation of ParticleSystem.
   * @param[in,out] ptc Particle system.
   * @param[in] step_size The step size of the integration.
   * @return true Collision detected.
   * @return false No collision detected.
   */
  template <typename ParticleSys>
  bool operator()(ParticleSys const &ptc, typename ParticleSys::Scalar step_size) {
    static_assert(HAS_METHOD(ParticleSys, radius), "Particle system doesn't have method radius()");

    auto const &px = ptc.pos().x;
    auto const &py = ptc.pos().y;
    auto const &pz = ptc.pos().z;
    auto const &r = ptc.radius();

    size_t num = ptc.number();
    for (size_t i = 0; i < num; ++i) {
      for (size_t j = i + 1; j < num; ++j) {
        auto dx = px[i] - px[j];
        auto dy = py[i] - py[j];
        auto dz = pz[i] - pz[j];

        if (sqrt(dx * dx + dy * dy + dz * dz) <= r[i] + r[j]) {
          return true;
        }
      }
    }
    return false;
  }
};

}  // namespace run_operations
#endif  // SPACEHUB_COLLISION_HPP
