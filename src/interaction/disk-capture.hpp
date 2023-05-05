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
 * @file diskcapture.hpp
 *
 * Header file.
 */
#pragma once

#include "../dev-tools.hpp"
#include "../spacehub-concepts.hpp"
namespace hub::force {
    class DiskCapture {
       public:
        constexpr static bool vel_dependent{true};
        // Type members
        template <typename Particles>
        static void add_acc_to(Particles const &particles, typename Particles::VectorArray &acceleration);

       private:
        template <typename Vec>
        double disk_rho(Vec &r, double M) {
            const double Q = 1;
            const double alpha = 1;
            const double lambda = 1;
            auto M_dot = lambda * 0.22 * (M / 1e7) / (2 * consts::pi);
            auto R = sqrt(r.x * r.x + r.y * r.y);
            auto Omega = sqrt(M / (R * R * R));
            auto H = pow(Q / 2 / alpha * M_dot / M / Omega, 1.0 / 3) * R;
            auto rho0 = M / (2 * Q * consts::pi * R * R * R);
            return rho0 * exp(-(r.z * r.z) / (2 * H * H));
        }

        template <typename Vec>
        Vec disk_v(Vec &r, double M) {
            auto rr = sqrt(r.x * r.x + r.y * r.y);
            auto v = sqrt(M / rr);
            return Vec{-v * r.y / rr, v * r.x / rr, 0};
        }
    };

    template <typename Particles>
    void DiskCapture::add_acc_to(const Particles &particles, typename Particles::VectorArray &acceleration) {
        size_t num = particles.number();
        auto const &p = particles.pos();
        auto const &v = particles.vel();
        auto const &m = particles.mass();
        auto const &r = particles.radius();

        for (size_t i = 1; i < num; ++i) {
            auto dr = p[i] - p[0];
            auto dv = v[i] - v[0];
            auto v_disk = disk_v(dr, m[0]);
            auto v_rel = dv - v_disk;
            auto rho = disk_rho(dr, m[0]);
            auto rd = r[i];
            auto vabs = sqrt(dot(v_rel, v_rel));
            /*double f1 = consts::pi * rd * rd * rho * dot(v_rel, v_rel);
            double f2 = 4 * consts::pi * consts::G * consts::G * m[i] * m[i] / dot(v_rel, v_rel) * rho;
            acceleration[i] -= std::max(f1, f2) * v_rel / vabs;
            acceleration[0] -= std::max(f1, f2) * v_rel / vabs;*/
            // double f1 = consts::pi * rd * rd * rho * dot(v_rel, v_rel);
            double f = 4 * consts::pi * consts::G * consts::G * m[i] * m[i] / dot(v_rel, v_rel) * rho;
            acceleration[i] -= f * v_rel / vabs / m[i];
            acceleration[0] += f * v_rel / vabs / m[0];
        }
    }

}  // namespace hub::force
