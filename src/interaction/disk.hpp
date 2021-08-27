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
 * @file disk.hpp
 *
 * Header file.
 */
#pragma once

#include "../dev-tools.hpp"
#include "../orbits/orbits.hpp"
#include "../spacehub-concepts.hpp"

namespace hub::force {

    class Disk {
       public:
        size_t const attached_BH_id{0};
        double const surface_density{0};
        double const R_in{0};
        double const R_out{0};
        double const h{0};
        double const T{0};
        double const opacity{0};
        vec3d const n{0, 0, 1};

        template <typename Vec3>
        bool is_in_disk(Vec3 const &dr) {
            // TBD
            return true;
        }
    };

    class DiskDamp {
       public:
        constexpr static bool vel_dependent{true};

        // Type members
        template <typename Particles>
        static void add_acc_to(Particles const &particles, typename Particles::VectorArray &acceleration);

        static void add_disk(Disk const &disk);

       private:
        template <typename Vector, typename Scalar>
        auto calc_t_damp(Disk const &disk, Scalar M, Scalar m, Vector const &dr, Vector const &dv);

        CREATE_METHOD_CHECK(chain_pos);

        CREATE_METHOD_CHECK(chain_vel);

        CREATE_METHOD_CHECK(index);

        static std::vector<Disk> disks;
    };

    void DiskDamp::add_disk(Disk const &disk) { disks.emplace_back(disk); }

    template <typename Vector, typename Scalar>
    auto DiskDamp::calc_t_damp(Disk const &disk, Scalar M, Scalar m, Vector const &dr, Vector const &dv) {
        auto [a, e, L] = orbit::calc_a_e_L(consts::G * (M + m), dr, dv);
        auto inc = acos(dot(L, disk.n) / (norm(L) * norm(disk.n)));
        // t_e t_i TBD
        return std::make_tuple(t_e, t_i);
    }

    template <typename Particles>
    void DiskDamp::add_acc_to(const Particles &particles, typename Particles::VectorArray &acceleration) {
        size_t num = particles.number();
        auto const &p = particles.pos();
        auto const &v = particles.vel();
        auto const &m = particles.mass();

        auto force = [&](auto const &dr, auto const &dv, auto i, auto j) {
            for (auto const &disk : disks) {
                if (i == disk.attached_BH_id && disk.is_in_disk(dr)) {
                    auto [t_e, t_i] = calc_t_damp(disk, m[i], m[j], dr, dv);

                    acceleration[j] += (-2 * dot(dv, dr) / norm2(dr) / t_e) * dr;  // damping in the radial direction
                    acceleration[j] += (-dot(dv, disk.n) / t_i) * disk.n;          // damping in the z of disk direction

                } else if (j == disk.attached_BH_id && disk.is_in_disk(dr)) {
                    auto [t_e, t_i] = calc_t_damp(disk, m[j], m[i], dr, dv);

                    acceleration[i] -= (-2 * dot(dv, dr) / norm2(dr) / t_e) * dr;
                    acceleration[i] -= (-dot(dv, disk.n) / t_i) * disk.n;
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

    class DiskTorque {
       public:
        constexpr static bool vel_dependent{true};

        // Type members
        template <typename Particles>
        static void add_acc_to(Particles const &particles, typename Particles::VectorArray &acceleration);

        static void add_disk(Disk const &disk);

       private:
       
        CREATE_METHOD_CHECK(chain_pos);

        CREATE_METHOD_CHECK(chain_vel);

        CREATE_METHOD_CHECK(index);

        static std::vector<Disk> disks;
    };

    void DiskTorque::add_disk(Disk const &disk) { disks.emplace_back(disk); }

    template <typename Particles>
    void DiskTorque::add_acc_to(const Particles &particles, typename Particles::VectorArray &acceleration) {
        size_t num = particles.number();
        auto const &p = particles.pos();
        auto const &v = particles.vel();
        auto const &m = particles.mass();

        auto force = [&](auto const &dr, auto const &dv, auto i, auto j) {
            for (auto const &disk : disks) {
                if (i == disk.attached_BH_id && disk.is_in_disk(dr)) {
                } else if (j == disk.attached_BH_id && disk.is_in_disk(dr)) {
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
}  // namespace hub::force
