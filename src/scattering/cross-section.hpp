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
 * @file cross-section.hpp
 *
 * Header file.
 */
#ifndef SPACEHUB_CROSS_SECTION_HPP
#define SPACEHUB_CROSS_SECTION_HPP

#include "../orbits/orbits.hpp"
#include "../rand-generator.hpp"
#include "../vector/vector3.hpp"

/**
 * @namespace space::scattering
 * namespace for scattering
 */
namespace space::scattering {

  template<typename Scalar>
  auto random_incident(Scalar m_stay, Scalar m_incident, Scalar v_inf, Scalar b_max, Scalar r) {
    using Vector = Vec3<Scalar>;
    auto b = sqrt(random::Uniform::get(0, b_max * b_max));
    auto w = random::Uniform::get(0, 2 * consts::pi);
    auto orbit = orbit::HyperOrbit(m_stay, m_incident, v_inf, b, r, w, 0, 0);
    Vector pos, vel;
    orbit::oribt_args_to_coord(orbit, pos, vel);
    return std::make_tuple(pos, vel);
  }

}

#endif //SPACEHUB_CROSS_SECTION_HPP
