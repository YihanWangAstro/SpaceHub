//
// Created by yihan on 2/25/19.
//

#ifndef SPACEHUB_NEWTONIAN_H
#define SPACEHUB_NEWTONIAN_H

#include "../interaction.h"

namespace SpaceH{
    template <typename TypeSystem>
    class NewtonianGrav : Interactions<NewtonianGrav<TypeSystem>, false>{
    public:
        SPACEHUB_USING_TYPE_SYSTEM_OF(TypeSystem);

        template<typename Particles>
        void impl_eval_acc(Particles const &partc, Coord &acc) {
            size_t num = partc.number();
            auto const &px = partc.pos().x;
            auto const &py = partc.pos().y;
            auto const &pz = partc.pos().z;
            auto const &m = partc.mass();

            calc::set_arrays_zero(acc.x, acc.y, acc.z);

            for (size_t i = 0; i < num; ++i) {
                for (size_t j = i + 1; j < num; ++j) {
                    auto dx = px[j] - px[i];
                    auto dy = py[j] - py[i];
                    auto dz = pz[j] - pz[i];
                    auto r = sqrt(dx * dx + dy * dy + dz * dz);
                    auto rr3 = 1.0 / (r * r * r);
                    acc.x[i] += dx * rr3 * m[j];
                    acc.y[i] += dy * rr3 * m[j];
                    acc.z[i] += dz * rr3 * m[j];
                    acc.x[j] -= dx * rr3 * m[i];
                    acc.y[j] -= dy * rr3 * m[i];
                    acc.z[j] -= dz * rr3 * m[i];
                }
            }
        }
    };
}
#endif //SPACEHUB_NEWTONIAN_H
