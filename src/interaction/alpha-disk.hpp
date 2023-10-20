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

#include <cmath>

#include "../dev-tools.hpp"
#include "../spacehub-concepts.hpp"
namespace hub::force {
    class AlphaDisk {
       public:
        constexpr static bool vel_dependent{true};
        // Type members
        template <typename Particles>
        static void add_acc_to(Particles const &particles, typename Particles::VectorArray &acceleration);

        static double Q;
        static double alpha;
        static double lambda;
        static int force;

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
    void AlphaDisk::add_acc_to(const Particles &particles, typename Particles::VectorArray &acceleration) {
        size_t num = particles.number();
        auto const &p = particles.pos();
        auto const &v = particles.vel();
        auto const &m = particles.mass();
        auto const &r = particles.radius();

        auto I_sup = [](double M, double logR) { return (0.5 * log(1 - 1 / M / M) + logR) / M / M; };
        auto I_sub = [](double M) { return (0.5 * log((1 + M) / (1 - M)) - M) / M / M; };
        auto dIdM_sup = [](double M, double logR) {
            return (-2 * logR + 1 / (M * M - 1) - log(1 - 1 / M / M)) / M / M / M;
        };
        auto dIdM_sub = [](double M) {
            return (M * M * M + (1 - M * M) * log((1 + M) / (1 - M)) - 2 * M) / (M * M * M * (M * M - 1));
        };

        constexpr double logR = 3.0;
        const double eps = 1 / std::exp(2.0 * logR / 3.0);
        const double x1 = 1 - eps;
        const double x2 = 1 + eps;
        const double y1 = I_sub(x1);
        const double y2 = I_sup(x2, logR);
        const double k1 = dIdM_sub(x1);
        const double k2 = dIdM_sup(x2, logR);
        const double a = k1 * (x2 - x1) - (y2 - y1);
        const double b = -k2 * (x2 - x1) + (y2 - y1);

        auto tt = [&](double M) { return (M - x1) / (x2 - x1); };

        auto connect = [&](double M) {
            double t = tt(M);
            return (1 - t) * y1 + y2 * t + (1 - t) * t * (t * b + (1 - t) * a);
        };

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

            double I = 0;

            double f = 0;

            if (force >= 0) {
                if (Mach > 1 + eps) {
                    I = (0.5 * log(1 - 1 / (Mach * Mach)) + logR) / (Mach * Mach);
                } else if ((0.1 < Mach) && (Mach < 1 - eps)) {
                    I = (0.5 * log((1 + Mach) / (1 - Mach)) - Mach) / (Mach * Mach);
                } else if (Mach < 0.1) {
                    I = Mach / 3.0;
                } else {
                    I = connect(Mach);
                }
            }

            if (force == 1) {
                f = I * 4 * consts::pi * consts::G * consts::G * m[i] * m[i] / (cs * cs) * rho;
            } else if (force == -1) {
                f = 4 * consts::pi * rd * rd * rho * v2;
            } else {
                double df = I * 4 * consts::pi * consts::G * consts::G * m[i] * m[i] / (cs * cs) * rho;
                double aeo = 4 * consts::pi * rd * rd * rho * v2;
                f = df + aeo;
            }
            acceleration[i] -= f * v_rel / vabs / m[i];
            acceleration[0] += f * v_rel / vabs / m[0];
        }
    }

}  // namespace hub::force
