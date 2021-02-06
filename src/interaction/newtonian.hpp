//
// Created by yihan on 2/25/19.
//

#pragma once

#include "../spacehub-concepts.hpp"
/**
 * @namespace space::interactions
 * Documentation for space
 */
namespace space::interactions {
    /*---------------------------------------------------------------------------*\
         Class NewtonianGrav Declaration
    \*---------------------------------------------------------------------------*/
    /**
     *
     */
    class NewtonianGrav {
       public:
        constexpr static bool vel_dependent{false};

        // Type members
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
        size_t num = particles.number();
        auto const &p = particles.pos();
        auto const &m = particles.mass();

        auto force = [&](auto const &dr, auto i, auto j) {
            auto r = norm(dr);
            auto rr3 = 1.0 / (r * r * r);
            auto mri = rr3 * m[i];
            auto mrj = rr3 * m[j];
            /*acceleration[i] += dr * mrj;
            acceleration[j] -= dr * mri;*/
            acceleration[i].x += dr.x * mrj;
            acceleration[i].y += dr.y * mrj;
            acceleration[i].z += dr.z * mrj;
            acceleration[j].x -= dr.x * mri;
            acceleration[j].y -= dr.y * mri;
            acceleration[j].z -= dr.z * mri;
        };

        if constexpr (HAS_METHOD(Particles, chain_pos) && HAS_METHOD(Particles, index)) {
            auto const &ch_p = particles.chain_pos();
            auto const &idx = particles.index();

            size_t size = particles.number();
            for (size_t i = 0; i < size - 1; ++i) {
                force(ch_p[i], idx[i], idx[i + 1]);
            }
            for (size_t i = 0; i < size - 2; ++i) {
                force(ch_p[i] + ch_p[i + 1], idx[i], idx[i + 2]);
            }

            for (size_t i = 0; i < size; ++i) {
                for (size_t j = i + 3; j < size; ++j) {
                    force(p[idx[j]] - p[idx[i]], idx[i], idx[j]);
                }
            }
        } else {
            for (size_t i = 0; i < num; ++i) {
                for (size_t j = i + 1; j < num; ++j) {
                    force(p[j] - p[i], i, j);
                }
            }
        }
    }
}  // namespace space::interactions
