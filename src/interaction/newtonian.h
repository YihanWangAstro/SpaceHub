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

        using Base = Interactions<NewtonianGrav<TypeSystem>, false>;
        using Base::ax;
        using Base::ay;
        using Base::az;

        template<typename Particles>
        void impl_eval_acc(Particles const &partc) {
            size_t num = partc.number();
            for (size_t i = 0; i < num; ++i) {
                for (size_t j = i + 1; j < num; ++j) {
                    auto dx  = partc.px[j] - partc.px[i];
                    auto dy  = partc.py[j] - partc.py[i];
                    auto dz  = partc.pz[j] - partc.pz[i];
                    auto r   = sqrt(dx * dx + dy * dy + dz * dz);
                    auto rr3 = 1 / (r * r * r);
                    ax[i] += dx * rr3 * partc.mass[j];
                    ay[i] += dy * rr3 * partc.mass[j];
                    az[i] += dz * rr3 * partc.mass[j];
                    ax[j] -= dx * rr3 * partc.mass[i];
                    ay[j] -= dy * rr3 * partc.mass[i];
                    az[j] -= dz * rr3 * partc.mass[i];
                }
            }
        }
    };
}
#endif //SPACEHUB_NEWTONIAN_H
