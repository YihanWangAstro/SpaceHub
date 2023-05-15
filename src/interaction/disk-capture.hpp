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

        static double Q;
        static double alpha;
        static double lambda;

        template <typename Vec>
        static auto local_Omega_H(Vec &r, double M) {
            auto M_dot = lambda * 0.22 * (M / 1e7) / (2 * consts::pi);
            auto R = sqrt(r.x * r.x + r.y * r.y);
            auto Omega = sqrt(M / (R * R * R));
            auto H = pow(Q / 2 / alpha * M_dot / M / Omega, 1.0 / 3) * R;
            return std::tuple(Omega, H);
        }

        template <typename Vec>
        static double disk_rho(Vec &r, double M, double H) {
            auto R = sqrt(r.x * r.x + r.y * r.y);
            auto rho0 = M / (2 * Q * consts::pi * R * R * R);
            return rho0 * exp(-(r.z * r.z) / (2 * H * H));
        }

        template <typename Vec>
        static Vec disk_v(Vec &r, double M) {
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

        constexpr double eps = 1e-3;
        constexpr double offset = log(2 * eps) / (1 + eps) / (1 + eps) * 0.5;
        // constexpr double logx = 0.5 * log((2 - eps) / eps * (1 + eps) * (1 + eps) / ((1 + eps) * (1 + eps) - 1)) - 1;

        for (size_t i = 1; i < num; ++i) {
            auto dr = p[i] - p[0];
            auto dv = v[i] - v[0];
            auto [Omega, H] = local_Omega_H(dr, m[0]);
            auto v_disk = disk_v(dr, m[0]);
            auto v_rel = dv - v_disk;
            auto rho = disk_rho(dr, m[0], H);
            auto rd = r[i];
            auto v2 = dot(v_rel, v_rel);
            auto vabs = sqrt(v2);
            auto cs = Omega * H;
            auto Mach = vabs / cs;
            auto logR = log(H / rd);

            double I = logR + offset;
            if (Mach > 1 + eps) {
                I = (0.5 * log(1 - 1 / (Mach * Mach)) + logR) / (Mach * Mach);
            } else if (Mach < 1 - eps) {
                I = (0.5 * log((1 + Mach) / (1 - Mach)) - Mach) / (Mach * Mach);
            }
            double df = I * 4 * consts::pi * consts::G * consts::G * m[i] * m[i] / (cs * cs) * rho;
            double aero_drag = 4 * consts::pi * rd * rd * rho * v2;
            double f = df + aero_drag;
            acceleration[i] -= f * v_rel / vabs / m[i];
            acceleration[0] += f * v_rel / vabs / m[0];
        }
    }

}  // namespace hub::force
