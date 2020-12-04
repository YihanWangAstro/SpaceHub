//
// Created by 王艺涵 on 4/12/19.
//

#pragma once

#include "../dev-tools.hpp"
#include "interaction.hpp"

namespace space::interactions {

    // pair-wise Post-Newtonian term in Centre of mass reference frame. on page 88.  https://arxiv.org/pdf/1310.1528.pdf
    class PN1 {
       public:
        // Type members
        template <typename Particles>
        static void add_acc_to(Particles const &particles, typename Particles::VectorArray &acceleration);

       private:
        CREATE_METHOD_CHECK(chain_pos);

        CREATE_METHOD_CHECK(index);
    };

    class PN2 {
       public:
        // Type members
        template <typename Particles>
        static void add_acc_to(Particles const &particles, typename Particles::VectorArray &acceleration);

       private:
        CREATE_METHOD_CHECK(chain_pos);

        CREATE_METHOD_CHECK(index);
    };

    class PN2p5 {
       public:
        // Type members
        template <typename Particles>
        static void add_acc_to(Particles const &particles, typename Particles::VectorArray &acceleration);

       private:
        CREATE_METHOD_CHECK(chain_pos);

        CREATE_METHOD_CHECK(index);
    };

    constexpr double INV_C = 1.0 / consts::C;
    constexpr double INV_C2 = INV_C * INV_C;
    constexpr double INV_C3 = INV_C * INV_C2;
    constexpr double INV_C4 = INV_C2 * INV_C2;
    constexpr double INV_C5 = INV_C2 * INV_C3;

    /*---------------------------------------------------------------------------*\
          Class PN1 Implememtation
    \*---------------------------------------------------------------------------*/
    template <typename Particles>
    void PN1::add_acc_to(const Particles &particles, typename Particles::VectorArray &acceleration) {
        size_t num = particles.number();
        auto &p = particles.pos();
        auto &v = particles.vel();
        auto &m = particles.mass();

        auto force = [&](auto const &dr, auto const &dv, auto i, auto j) {
            auto m_tot = m[i] + m[j];
            auto m_nu = m[i] * m[j] / (m_tot * m_tot);

            auto r2 = norm2(dr);
            auto r = sqrt(r2);
            auto n = dr / r;
            auto r_dot = dot(n, dr);

            auto v2 = norm2(dv);

            auto A1 = -1.5 * r_dot * r_dot * m_nu + v2 + 3 * m_nu * v2;
            auto B = r_dot * (-4 + 2 * m_nu);
            auto A2 = -consts::G / r * (4 + 2 * m_nu);

            auto Ai = A1 + m[j] * A2;
            auto Aj = A1 + m[i] * A2;

            auto coef = consts::G / r2 * INV_C2;

            acceleration[i] += (coef * m[j]) * (Ai * n + B * dv);
            acceleration[j] -= (coef * m[i]) * (Aj * n + B * dv);
        };

        if constexpr (HAS_METHOD(Particles, chain_pos) && HAS_METHOD(Particles, index)) {
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
              Class PN2 Implememtation
    \*---------------------------------------------------------------------------*/
    template <typename Particles>
    void PN2::add_acc_to(const Particles &particles, typename Particles::VectorArray &acceleration) {
        size_t num = particles.number();
        auto &p = particles.pos();
        auto &v = particles.vel();
        auto &m = particles.mass();

        auto force = [&](auto const &dr, auto const &dv, auto i, auto j) {
            auto m_tot = m[i] + m[j];

            auto m_nu = m[i] * m[j] / (m_tot * m_tot);
            auto m_nu2 = m_nu * m_nu;

            auto r2 = norm2(dr);
            auto r = sqrt(r2);
            auto n = dr / r;
            auto r_dot = dot(n, dr);
            auto r_dot2 = r_dot * r_dot;
            auto r_dot3 = r_dot2 * r_dot;
            auto r_dot4 = r_dot2 * r_dot2;

            auto v2 = norm2(dv);

            auto gmri = consts::G * m[i] / r;
            auto gmri2 = gmri * gmri;
            auto gmrj = consts::G * m[j] / r;
            auto gmrj2 = gmrj * gmrj;

            auto A1 = r_dot4 * (1.875 * m_nu - 5.625 * m_nu2) - r_dot2 * (4.5 * m_nu * v2 + 6 * m_nu2 * v2) +
                      v2 * v2 * (3 * m_nu - 4 * m_nu2);

            auto A2 = r_dot2 * (-2 - 25 * m_nu - 2 * m_nu2) - v2 * (6.5 * m_nu + 2 * m_nu2);

            auto A3 = 9 + 21.75 * m_nu;

            auto B1 = 4.5 * r_dot3 * m_nu + 3 * r_dot3 * m_nu2 - 7.5 * r_dot * m_nu * v2 - 2 * r_dot * m_nu2 * v2;

            auto B2 = r_dot * (2 + 20.5 * m_nu + 4 * m_nu2);

            auto Ai = A1 + gmrj * A2 + gmrj2 * A3;

            auto Aj = A1 * gmri * A2 + gmri2 * A3;

            auto Bi = B1 + B2 * gmrj;

            auto Bj = B1 + B2 * gmri;

            auto coef = consts::G / r2 * INV_C4;

            acceleration[i] += (coef * m[j]) * (Ai * n + Bi * dv);
            acceleration[j] -= (coef * m[i]) * (Aj * n + Bj * dv);
        };

        if constexpr (HAS_METHOD(Particles, chain_pos) && HAS_METHOD(Particles, index)) {
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
          Class PN2.5 Implememtation
    \*---------------------------------------------------------------------------*/
    template <typename Particles>
    void PN2p5::add_acc_to(const Particles &particles, typename Particles::VectorArray &acceleration) {
        size_t num = particles.number();
        auto &p = particles.pos();
        auto &v = particles.vel();
        auto &m = particles.mass();

        auto force = [&](auto const &dr, auto const &dv, auto i, auto j) {
            auto m_tot = m[i] + m[j];
            auto m_nu = m[i] * m[j] / (m_tot * m_tot);

            auto r2 = norm2(dr);
            auto r = sqrt(r2);
            auto n = dr / r;
            auto r_dot = dot(n, dr);

            auto v2 = norm2(dv);

            auto gmri = consts::G * m[i] / r;
            auto gmrj = consts::G * m[j] / r;

            auto Ai = gmrj * m_nu * r_dot * (-4.8 * v2 - 136.0 / 15 * gmrj);
            auto Aj = gmri * m_nu * r_dot * (-4.8 * v2 - 136.0 / 15 * gmri);

            auto Bi = gmrj * m_nu * (1.6 * v2 + 4.8 * gmrj);
            auto Bj = gmri * m_nu * (1.6 * v2 + 4.8 * gmri);

            auto coef = consts::G / r2 * INV_C5;

            acceleration[i] += (coef * m[j]) * (Ai * n + Bi * dv);
            acceleration[j] -= (coef * m[i]) * (Aj * n + Bj * dv);
        };

        if constexpr (HAS_METHOD(Particles, chain_pos) && HAS_METHOD(Particles, index)) {
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
}  // namespace space::interactions
