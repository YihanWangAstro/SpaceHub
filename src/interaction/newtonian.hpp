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
 * @file newtonian.hpp
 *
 * Header file.
 */
#pragma once

#include "../spacehub-concepts.hpp"
/**
 * @namespace space::force
 * Documentation for space
 */
namespace space::force {
    /*---------------------------------------------------------------------------*\
         Class NewtonianGrav Declaration
    \*---------------------------------------------------------------------------*/

    /**
     * @brief Newtonian direct summation force.
     *
     */
    class NewtonianGrav {
       public:
        /**
         * @brief Is this force velocity dependent?
         *
         */
        constexpr static bool vel_dependent{false};

        // Type members
        /**
         * @brief Add newtonian acceleration to existing 3D vector array.
         *
         * @note this method ADD newtonian acceleration TO input 'acceleration'.
         *
         * @tparam Particles Particle system type satisfy concept particle system.
         * @param[in] particles Particle system that is used to evaluated the acceleration.
         * @param[in,out] acceleration 3D vector array to be updated.
         */
        template <typename Particles>
        static void add_acc_to(Particles const &particles, typename Particles::VectorArray &acceleration);

       private:
        CREATE_METHOD_CHECK(chain_pos);

        CREATE_METHOD_CHECK(index);
    };

    /*---------------------------------------------------------------------------*\
          Class NewtonianGrav Implementation
    \*---------------------------------------------------------------------------*/
    template <typename Particles>
    void NewtonianGrav::add_acc_to(const Particles &particles, typename Particles::VectorArray &acceleration) {
        using Vector = typename Particles::Vector;
        size_t num = particles.number();
        auto const &p = particles.pos();
        auto const &m = particles.mass();

        auto force = [&](Vector const &dr, size_t i, size_t j) {
            auto r = norm(dr);
            auto rr3 = 1.0 / (r * r * r);
            /*
            acceleration[i] += dr * rr3 * m[j];
            acceleration[j] -= dr * rr3 * m[i];*/
            acceleration[i].x += dr.x * rr3 * m[j];
            acceleration[i].y += dr.y * rr3 * m[j];
            acceleration[i].z += dr.z * rr3 * m[j];
            acceleration[j].x -= dr.x * rr3 * m[i];
            acceleration[j].y -= dr.y * rr3 * m[i];
            acceleration[j].z -= dr.z * rr3 * m[i];
        };

        if constexpr (HAS_METHOD(Particles, chain_pos) && HAS_METHOD(Particles, index)) {
            auto const &ch_p = particles.chain_pos();
            auto const &idx = particles.index();

            size_t size = particles.number();

            for (size_t i = 0; i < size; ++i) {
                for (size_t j = i + 3; j < size; ++j) {
                    force(p[idx[j]] - p[idx[i]], idx[i], idx[j]);
                }
            }

            for (size_t i = 0; i < size - 2; ++i) {
                force(ch_p[i] + ch_p[i + 1], idx[i], idx[i + 2]);
            }

            for (size_t i = 0; i < size - 1; ++i) {
                force(ch_p[i], idx[i], idx[i + 1]);
            }
        } else {
            for (size_t i = 0; i < num; ++i) {
                for (size_t j = i + 1; j < num; ++j) {
                    force(p[j] - p[i], i, j);
                }
            }
        }
    }
}  // namespace space::force
