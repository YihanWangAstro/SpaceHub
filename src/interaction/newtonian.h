//
// Created by yihan on 2/25/19.
//

#ifndef SPACEHUB_NEWTONIAN_H
#define SPACEHUB_NEWTONIAN_H

#include "../interaction.h"
#include "../dev-tools.h"

namespace space::interactions {
    /*---------------------------------------------------------------------------*\
        Class NewtonianGrav Declaration
    \*---------------------------------------------------------------------------*/
    class NewtonianGrav : public Interactions<NewtonianGrav> {
    public:
        //Public methods
        template<typename Particles>
        void impl_eval_acc(Particles const &partc, typename Particles::Coord &acc);

        template<typename Particles>
        void impl_eval_newtonian_acc(Particles const &partc, typename Particles::Coord &acc);

    private:
        CREATE_METHOD_CHECK(chain_pos);

        CREATE_METHOD_CHECK(chain_vel);

        CREATE_METHOD_CHECK(index);
    };

    /*---------------------------------------------------------------------------*\
        Class NewtonianGrav Implememtation
    \*---------------------------------------------------------------------------*/
    template<typename Particles>
    void NewtonianGrav::impl_eval_acc(const Particles &partc, typename Particles::Coord &acc) {
        impl_eval_newtonian_acc(partc, acc);
    }

    template<typename Particles>
    void NewtonianGrav::impl_eval_newtonian_acc(const Particles &partc, typename Particles::Coord &acc) {
        size_t num = partc.number();
        auto &px = partc.pos().x;
        auto &py = partc.pos().y;
        auto &pz = partc.pos().z;
        auto &m = partc.mass();

        calc::set_arrays_zero(acc.x, acc.y, acc.z);

        auto force = [&](auto dx, auto dy, auto dz, auto i, auto j) {
            auto r = sqrt(dx * dx + dy * dy + dz * dz);
            auto rr3 = 1.0 / (r * r * r);
            acc.x[i] += dx * rr3 * m[j];
            acc.y[i] += dy * rr3 * m[j];
            acc.z[i] += dz * rr3 * m[j];
            acc.x[j] -= dx * rr3 * m[i];
            acc.y[j] -= dy * rr3 * m[i];
            acc.z[j] -= dz * rr3 * m[i];
        };

        if constexpr (HAS_METHOD(Particles, chain_pos) && HAS_METHOD(Particles, chain_vel) &&
                      HAS_METHOD(Particles, index)) {
            auto const &ch_px = partc.chain_pos().x;
            auto const &ch_py = partc.chain_pos().y;
            auto const &ch_pz = partc.chain_pos().z;
            auto const &ch_vx = partc.chain_vel().x;
            auto const &ch_vy = partc.chain_vel().y;
            auto const &ch_vz = partc.chain_vel().z;
            auto const &idx = partc.index();

            size_t size = partc.number();
            for (size_t i = 0; i < size - 1; ++i)
                force(ch_px[i], ch_py[i], ch_pz[i], idx[i], idx[i + 1]);

            for (size_t i = 0; i < size - 2; ++i)
                force(ch_px[i] + ch_px[i + 1], ch_py[i] + ch_py[i + 1], ch_pz[i] + ch_pz[i + 1], idx[i], idx[i + 2]);

            for (size_t i = 0; i < size; ++i)
                for (size_t j = i + 3; j < size; ++j)
                    force(px[idx[j]] - px[idx[i]], py[idx[j]] - py[idx[i]], pz[idx[j]] - pz[idx[i]], idx[i], idx[j]);

        } else {
            for (size_t i = 0; i < num; ++i) {
                for (size_t j = i + 1; j < num; ++j) {
                    force(px[j] - px[i], py[j] - py[i], pz[j] - pz[i], i, j);
                }
            }
        }

    }
}
#endif //SPACEHUB_NEWTONIAN_H
