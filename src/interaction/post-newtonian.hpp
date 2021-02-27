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
#include "interaction.hpp"

namespace space::force {

    // pair-wise Post-Newtonian term in Centre of mass reference frame. on page 76. https://arxiv.org/pdf/1310.1528.pdf

    /**
     * @brief First order Post-Newtonian term for general relativity.
     *
     */
    class PN1 {
       public:
        /**
         * @brief Is this force velocity dependent?
         *
         */
        constexpr static bool vel_dependent{true};

        // Type members
        /**
         * @brief Add acceleration from first order post-newtonian term to existing 3D vector array.
         *
         * @note this method ADD acceleration TO input 'acceleration'.
         *
         * @tparam Particles Particle system type satisfy concept particle system.
         * @param[in] particles Particle system that is used to evaluated the acceleration.
         * @param[in,out] acceleration 3D vector array to be updated.
         */
        template <typename Particles>
        static void add_acc_to(Particles const &particles, typename Particles::VectorArray &acceleration);

       private:
        CREATE_METHOD_CHECK(chain_pos);

        CREATE_METHOD_CHECK(chain_vel);

        CREATE_METHOD_CHECK(index);
    };

    class PN2 {
       public:
        /**
         * @brief Is this force velocity dependent?
         *
         */
        constexpr static bool vel_dependent{true};

        // Type members
        /**
         * @brief Add acceleration from second order post-newtonian term to existing 3D vector array.
         *
         * @note this method ADD acceleration TO input 'acceleration'.
         *
         * @tparam Particles Particle system type satisfy concept particle system.
         * @param[in] particles Particle system that is used to evaluated the acceleration.
         * @param[in,out] acceleration 3D vector array to be updated.
         */
        template <typename Particles>
        static void add_acc_to(Particles const &particles, typename Particles::VectorArray &acceleration);

       private:
        CREATE_METHOD_CHECK(chain_pos);

        CREATE_METHOD_CHECK(chain_vel);

        CREATE_METHOD_CHECK(index);
    };

    class PN2p5 {
       public:
        /**
         * @brief Is this force velocity dependent?
         *
         */
        constexpr static bool vel_dependent{true};

        // Type members
        /**
         * @brief Add acceleration from two point five order post-newtonian term to existing 3D vector array.
         *
         * @note this method ADD acceleration TO input 'acceleration'.
         *
         * @tparam Particles Particle system type satisfy concept particle system.
         * @param[in] particles Particle system that is used to evaluated the acceleration.
         * @param[in,out] acceleration 3D vector array to be updated.
         */
        template <typename Particles>
        static void add_acc_to(Particles const &particles, typename Particles::VectorArray &acceleration);

       private:
        CREATE_METHOD_CHECK(chain_pos);

        CREATE_METHOD_CHECK(chain_vel);

        CREATE_METHOD_CHECK(index);
    };

    constexpr double INV_C = 1.0 / consts::C;
    constexpr double INV_C2 = INV_C * INV_C;
    constexpr double INV_C3 = INV_C * INV_C2;
    constexpr double INV_C4 = INV_C2 * INV_C2;
    constexpr double INV_C5 = INV_C2 * INV_C3;

    /*---------------------------------------------------------------------------*\
          Class PN1 Implementation
    \*---------------------------------------------------------------------------*/
    template <typename Particles>
    void PN1::add_acc_to(const Particles &particles, typename Particles::VectorArray &acceleration) {
        size_t num = particles.number();
        auto const &p = particles.pos();
        auto const &v = particles.vel();
        auto const &m = particles.mass();

        auto force = [&](auto const &dr, auto const &dv, auto i, auto j) {
            auto r2 = norm2(dr);
            auto r = sqrt(r2);
            auto n = -dr / r;

            auto v1s = norm2(v[i]);
            auto v2s = norm2(v[j]);
            auto v12 = dot(v[i], v[j]);

            auto nv1 = dot(n, v[i]);
            auto nv2 = dot(n, v[j]);

            auto gmr1 = consts::G * m[i] / r;
            auto gmr2 = consts::G * m[j] / r;

            auto Ai = -v1s - 2 * v2s + 4 * v12 + 1.5 * nv2 * nv2 + 5 * gmr1 + 4 * gmr2;

            auto Aj = -v2s - 2 * v1s + 4 * v12 + 1.5 * nv1 * nv1 + 5 * gmr2 + 4 * gmr1;

            auto Bi = 4 * nv1 - 3 * nv2;

            auto Bj = -4 * nv2 + 3 * nv1;

            auto coef = consts::G / r2 * INV_C2;

            acceleration[i] += (coef * m[j]) * (Ai * n - Bi * dv);
            acceleration[j] -= (coef * m[i]) * (Aj * n - Bj * dv);
        };

        if constexpr (HAS_METHOD(Particles, chain_pos) && HAS_METHOD(Particles, index) &&
                      HAS_METHOD(Particles, chain_vel)) {
            auto const &ch_p = particles.chain_pos();
            auto const &ch_v = particles.chain_vel();
            auto const &idx = particles.index();

            size_t size = particles.number();
            for (size_t i = 0; i < size - 1; ++i) {
                force(ch_p[i], ch_v[i], idx[i], idx[i + 1]);
            }

            for (size_t i = 0; i < size - 2; ++i) {
                force(ch_p[i] + ch_p[i + 1], ch_v[i] + ch_v[i + 1], idx[i], idx[i + 2]);
            }

            for (size_t i = 0; i < size; ++i) {
                for (size_t j = i + 3; j < size; ++j) {
                    force(p[idx[j]] - p[idx[i]], v[idx[j]] - v[idx[i]], idx[i], idx[j]);
                }
            }

        } else {
            for (size_t i = 0; i < num; ++i) {
                for (size_t j = i + 1; j < num; ++j) {
                    force(p[j] - p[i], v[j] - v[i], i, j);
                }
            }
        }
    }

    /*---------------------------------------------------------------------------*\
              Class PN2 Implementation
    \*---------------------------------------------------------------------------*/
    template <typename Particles>
    void PN2::add_acc_to(const Particles &particles, typename Particles::VectorArray &acceleration) {
        size_t num = particles.number();
        auto const &p = particles.pos();
        auto const &v = particles.vel();
        auto const &m = particles.mass();

        auto force = [&](auto const &dr, auto const &dv, auto i, auto j) {
            auto r2 = norm2(dr);
            auto r = sqrt(r2);
            auto n = -dr / r;

            auto v1s = norm2(v[i]);
            auto v1q = v1s * v1s;
            auto v2s = norm2(v[j]);
            auto v2q = v2s * v2s;
            auto v12 = dot(v[i], v[j]);

            auto nv1 = dot(n, v[i]);
            auto nv1s = nv1 * nv1;
            auto nv2 = dot(n, v[j]);
            auto nv2s = nv2 * nv2;

            auto gmr1 = consts::G * m[i] / r;
            auto gmr2 = consts::G * m[j] / r;

            auto m1s = m[i] * m[i];
            auto m2s = m[j] * m[j];
            auto m12 = m[i] * m[j];

            auto Ai = -2 * v2q + 4 * v2s * v12 - 2 * v12 * v12 +
                      nv2s * (1.5 * v1s + 4.5 * v2s - 6 * v12 - 1.875 * nv2s) +
                      gmr1 * (-3.75 * v1s + 1.25 * v2s - 2.5 * v12 + 19.5 * nv1s - 39 * nv1 * nv2 + 8.5 * nv2s) +
                      gmr2 * (4 * v2s - 8 * v12 + 2 * nv1s - 4 * nv1 * nv2 - 6 * nv2s) +
                      consts::G * consts::G / r2 * (-14.25 * m1s - 9 * m2s - 34.5 * m12);

            auto Aj = -2 * v1q + 4 * v1s * v12 - 2 * v12 * v12 +
                      nv1s * (1.5 * v2s + 4.5 * v1s - 6 * v12 - 1.875 * nv1s) +
                      gmr2 * (-3.75 * v2s + 1.25 * v1s - 2.5 * v12 + 19.5 * nv2s - 39 * nv1 * nv2 + 8.5 * nv1s) +
                      gmr1 * (4 * v1s - 8 * v12 + 2 * nv2s - 4 * nv1 * nv2 - 6 * nv1s) +
                      consts::G * consts::G / r2 * (-14.25 * m2s - 9 * m1s - 34.5 * m12);

            auto Bi = v1s * nv2 + 4 * v2s * nv1 - 5 * v2s * nv2 - 4 * v12 * nv1 + 4 * v12 * nv2 - 6 * nv1 * nv2s +
                      4.5 * nv2 * nv2s + gmr1 * (-15.75 * nv1 + 13.75 * nv2) + gmr2 * (-2 * nv1 - 2 * nv2);

            auto Bj = -v2s * nv1 - 4 * v1s * nv2 + 5 * v1s * nv2 + 4 * v12 * nv2 - 4 * v12 * nv1 + 6 * nv2 * nv1s -
                      4.5 * nv1 * nv1s + gmr2 * (15.75 * nv2 - 13.75 * nv1) + gmr1 * (2 * nv2 + 2 * nv1);

            auto coef = consts::G / r2 * INV_C4;

            acceleration[i] += (coef * m[j]) * (Ai * n - Bi * dv);
            acceleration[j] -= (coef * m[i]) * (Aj * n - Bj * dv);
        };

        if constexpr (HAS_METHOD(Particles, chain_pos) && HAS_METHOD(Particles, index) &&
                      HAS_METHOD(Particles, chain_vel)) {
            auto const &ch_p = particles.chain_pos();
            auto const &ch_v = particles.chain_vel();
            auto const &idx = particles.index();

            size_t size = particles.number();
            for (size_t i = 0; i < size - 1; ++i) {
                force(ch_p[i], ch_v[i], idx[i], idx[i + 1]);
            }

            for (size_t i = 0; i < size - 2; ++i) {
                force(ch_p[i] + ch_p[i + 1], ch_v[i] + ch_v[i + 1], idx[i], idx[i + 2]);
            }

            for (size_t i = 0; i < size; ++i) {
                for (size_t j = i + 3; j < size; ++j) {
                    force(p[idx[j]] - p[idx[i]], v[idx[j]] - v[idx[i]], idx[i], idx[j]);
                }
            }
        } else {
            for (size_t i = 0; i < num; ++i) {
                for (size_t j = i + 1; j < num; ++j) {
                    force(p[j] - p[i], v[j] - v[i], i, j);
                }
            }
        }
    }

    /*---------------------------------------------------------------------------*\
          Class PN2.5 Implementation
    \*---------------------------------------------------------------------------*/
    template <typename Particles>
    void PN2p5::add_acc_to(const Particles &particles, typename Particles::VectorArray &acceleration) {
        size_t num = particles.number();
        auto const &p = particles.pos();
        auto const &v = particles.vel();
        auto const &m = particles.mass();

        auto force = [&](auto const &dr, auto const &dv, auto i, auto j) {
            auto r2 = norm2(dr);
            auto r = sqrt(r2);
            auto n = -dr / r;

            auto v1s = norm2(v[i]);
            auto v2s = norm2(v[j]);
            auto v12 = dot(v[i], v[j]);
            auto dv2 = norm2(dv);

            /*auto nv1 = dot(n, v[i]);
            auto nv2 = dot(n, v[j]);*/

            auto nv = -dot(n, dv);

            auto gmr1 = consts::G * m[i] / r;
            auto gmr2 = consts::G * m[j] / r;

            auto Ai = nv * (3 * dv2 - 6 * gmr1 + 52.0 / 3 * gmr2);

            auto Aj = nv * (3 * dv2 - 6 * gmr2 + 52.0 / 3 * gmr1);

            auto Bi = -dv2 + 2 * gmr1 - 8 * gmr2;

            auto Bj = -dv2 + 2 * gmr2 - 8 * gmr1;

            auto coef = 0.8 * consts::G * consts::G * m[i] * m[j] / (r2 * r) * INV_C5;

            acceleration[i] += coef * (Ai * n - Bi * dv);
            acceleration[j] -= coef * (Aj * n - Bj * dv);
        };

        if constexpr (HAS_METHOD(Particles, chain_pos) && HAS_METHOD(Particles, index) &&
                      HAS_METHOD(Particles, chain_vel)) {
            auto const &ch_p = particles.chain_pos();
            auto const &ch_v = particles.chain_vel();
            auto const &idx = particles.index();

            size_t size = particles.number();
            for (size_t i = 0; i < size - 1; ++i) {
                force(ch_p[i], ch_v[i], idx[i], idx[i + 1]);
            }

            for (size_t i = 0; i < size - 2; ++i) {
                force(ch_p[i] + ch_p[i + 1], ch_v[i] + ch_v[i + 1], idx[i], idx[i + 2]);
            }

            for (size_t i = 0; i < size; ++i) {
                for (size_t j = i + 3; j < size; ++j) {
                    force(p[idx[j]] - p[idx[i]], v[idx[j]] - v[idx[i]], idx[i], idx[j]);
                }
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
