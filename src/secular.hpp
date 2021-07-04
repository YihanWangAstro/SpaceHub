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
 * @file secular.hpp
 *
 * Header file.
 */
#pragma once
#include "macros.hpp"

namespace hub::secular {

    template <typename T = double>
    T GW_dadt() {}

    template <typename T = double>
    T radial_tidal_dadt(T m, T M, T R, T a, T e, T k, T tau) {}

    template <typename T = double>
    T ELK_quad_timescale(T m1, T m2, T m3, T a1, T a2, T e2) {
        T m_in = m1 + m2;
        T P = 2 * consts::pi * sqrt(a1 * a1 * a1 / consts::G * m_in);
        T a_ratio_eff = a2 * sqrt(1 - e2 * e2) / a1;
        return P * m_in / m3 * a_ratio_eff * a_ratio_eff * a_ratio_eff;
    }

    template <typename EllipticOrb>
    auto ELK_quad_timescale(EllipticOrb const& orb1, EllipticOrb const& orb2) {
        return ELK_quad_timescale(orb1.m1, orb1.m2, orb2.m2, orb1.a, orb2.a, orb2.e);
    }

}  // namespace hub::secular