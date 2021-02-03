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
 * @file macros.hpp
 *
 * Header file.
 */
#pragma once

/**
 * @namespace space::consts
 * Documentation for space
 */
namespace space::consts {
    constexpr double pi = 3.14159265358979323846;
}

/**
 * @namespace space::unit
 * Documentation for space
 */
namespace space::unit {

#define MAKE_UNIT(NAME, UNIT)                                                                                 \
    constexpr double NAME = UNIT;                                                                             \
    constexpr double operator"" _##NAME(unsigned long long int x) { return static_cast<double>(x) * (UNIT); } \
    constexpr double operator"" _##NAME(long double x) { return static_cast<double>(x) * (UNIT); }

    MAKE_UNIT(AU, 1)
    MAKE_UNIT(Ms, 1)
    MAKE_UNIT(year, 2 * consts::pi)
    MAKE_UNIT(deg, consts::pi / 180.0)

    MAKE_UNIT(kyr, 1e3_year)
    MAKE_UNIT(Myr, 1e6_year)
    MAKE_UNIT(Gyr, 1e9_year)
    MAKE_UNIT(month, 1_year / 12.0)
    MAKE_UNIT(day, 1_year / 365.25636042)
    MAKE_UNIT(hr, 1_day / 24.0)
    MAKE_UNIT(min, 1_hr / 60.0)
    MAKE_UNIT(sec, 1_min / 60.0)
    MAKE_UNIT(T_hubble, 13.7_Gyr)

    MAKE_UNIT(km, 1_AU / 149597870.7)
    MAKE_UNIT(PC, 1_AU * 648000.0 / consts::pi)
    MAKE_UNIT(Rs, 6.957e5_km)
    MAKE_UNIT(Rj, 69911_km)

    MAKE_UNIT(kms, km / sec)

    MAKE_UNIT(Me, 3.003E-6_Ms)
    MAKE_UNIT(Mj, 317.8_Me)
    MAKE_UNIT(Mmoon, 0.012300_Me)

}  // namespace space::unit

namespace space::consts {
    constexpr double a_jupiter = 5.2044 * unit::AU;
    constexpr double e_jupiter = 0.0489;
    constexpr double LoAN_jupiter = 100.464 * unit::deg;
    constexpr double AoP_jupiter = 273.867 * unit::deg;
    constexpr double i_jupiter = 6.09 * unit::deg;

    constexpr double a_neptune = 30.11 * unit::AU;
    constexpr double e_neptune = 0.009456;
    constexpr double LoAN_neptune = 131.784 * unit::deg;
    constexpr double AoP_neptune = 276.336 * unit::deg;
    constexpr double i_neptune = 6.43 * unit::deg;

    constexpr double a_mercury = 0.387098 * unit::AU;
    constexpr double e_mercury = 0.20563;
    constexpr double LoAN_mercury = 48.331 * unit::deg;
    constexpr double AoP_mercury = 29.124 * unit::deg;
    constexpr double i_mercury = 3.38 * unit::deg;
}  // namespace space::consts

namespace space::consts {
    constexpr double G = 1;
    constexpr double C = 299792.458 * unit::kms;
}  // namespace space::consts
