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
 * @file post-newtonian.hpp
 *
 * Header file.
 */
#pragma once

#include "../dev-tools.hpp"
#include "../spacehub-concepts.hpp"
namespace space::force {
    class Tidal {
       public:
        constexpr static bool vel_dependent{true};

        // Type members
        template <typename Particles>
        static void add_acc_to(Particles const &particles, typename Particles::VectorArray &acceleration);

       private:
        CREATE_METHOD_CHECK(chain_pos);

        CREATE_METHOD_CHECK(chain_vel);

        CREATE_METHOD_CHECK(index);
    };

    template <typename Particles>
    void Tidal::add_acc_to(const Particles &particles, typename Particles::VectorArray &acceleration) {
        size_t num = particles.number();
        auto const &p = particles.pos();
        auto const &v = particles.vel();
        auto const &m = particles.mass();
        auto const &k = particles.tide_apsidal_const();
        auto const &tau = particles.tide_lag_time();
        auto const &rad = particles.radius();

        auto force = [&](auto const &dr, auto const &dv, auto i, auto j) {
            if (k[i] != 0 || k[j] != 0) {
                auto r = norm(dr);
                auto r2 = r * r;
                auto r4 = r2 * r2;
                auto r8 = r4 * r4;
                if (k[i] != 0) {
                    auto rad2 = rad[i] * rad[i];
                    auto rad4 = rad2 * rad2;
                    auto rad5 = rad4 * rad[i];
                    auto coef = 3 * consts::G * k[i] * rad5 * m[j] * m[j] * (1 + 3 * tau[i] * dot(dv, dr) / r2) / r8;
                    acceleration[i] += coef / m[i] * dr;
                    acceleration[j] -= coef / m[j] * dr;
                }
                if (k[j] != 0) {
                    auto rad2 = rad[j] * rad[j];
                    auto rad4 = rad2 * rad2;
                    auto rad5 = rad4 * rad[j];
                    auto coef = 3 * consts::G * k[j] * rad5 * m[i] * m[i] * (1 + 3 * tau[j] * dot(dv, dr) / r2) / r8;
                    acceleration[i] += coef / m[i] * dr;
                    acceleration[j] -= coef / m[j] * dr;
                }
            }
        };

        if constexpr (HAS_METHOD(Particles, chain_pos) && HAS_METHOD(Particles, index) &&
                      HAS_METHOD(Particles, chain_vel)) {
            auto const &ch_p = particles.chain_pos();
            auto const &ch_v = particles.chain_vel();
            auto const &idx = particles.index();

            size_t size = particles.number();

            for (size_t i = 0; i < size; ++i) {
                for (size_t j = i + 3; j < size; ++j) {
                    force(p[idx[j]] - p[idx[i]], v[idx[j]] - v[idx[i]], idx[i], idx[j]);
                }
            }

            for (size_t i = 0; i < size - 2; ++i) {
                force(ch_p[i] + ch_p[i + 1], ch_v[i] + ch_v[i + 1], idx[i], idx[i + 2]);
            }

            for (size_t i = 0; i < size - 1; ++i) {
                force(ch_p[i], ch_v[i], idx[i], idx[i + 1]);
            }
        } else {
            for (size_t i = 0; i < num; ++i) {
                for (size_t j = i + 1; j < num; ++j) {
                    force(p[j] - p[i], v[j] - v[i], i, j);
                }
            }
        }
    }
}  // namespace space::force
