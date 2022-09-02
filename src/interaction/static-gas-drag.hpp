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
 * @file gas-drag.hpp
 *
 * Header file.
 */
#pragma once

#include "../dev-tools.hpp"
#include "../spacehub-concepts.hpp"
namespace hub::force {
    class StaticGasDrag {
       public:
        constexpr static bool vel_dependent{true};

        // Type members
        template <typename Particles>
        static void add_acc_to(Particles const &particles, typename Particles::VectorArray &acceleration);

        static double subsonic_coef;
        static double sound_speed;
    };

    template <typename Particles>
    void StaticGasDrag::add_acc_to(const Particles &particles, typename Particles::VectorArray &acceleration) {
        size_t num = particles.number();
        auto const &v = particles.vel();
        auto const &m = particles.mass();
        const double supersonic_coef =
            3 * sqrt(consts::pi / 2) * subsonic_coef * sound_speed * sound_speed * sound_speed;

        for (size_t i = 0; i < num; ++i) {
            auto u = norm(v[i]);
            if (u < sound_speed) {
                acceleration[i] -= (subsonic_coef * m[i]) * v[i];
            } else {
                acceleration[i] -= (subsonic_coef * m[i] / (u * u * u)) * v[i];
            }
        }
    }

}  // namespace hub::force
