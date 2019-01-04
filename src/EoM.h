//
// Created by yihan on 1/3/19.
//

#ifndef SPACEHUB_EOM_H
#define SPACEHUB_EOM_H

#include "type_class.h"

namespace SpaceH{

    /*size_t num = partc.number();
    for(size_t i = 0 ; i < num; ++i){
    for(size_t j = i+1 ; j < num; ++j){
    auto dx  = partc.px(j) - partc.px(i);
    auto dy  = partc.py(j) - partc.py(i);
    auto dz  = partc.pz(j) - partc.pz(i);
    auto r   = sqrt(dx*dx+dy*dy+dz*dz);
    auto rr3 = 1/(r*r*r);
    acc.ax(i) +=  dx*rr3*partc.mass(j);
    acc.ay(i) +=  dy*rr3*partc.mass(j);
    acc.az(i) +=  dz*rr3*partc.mass(j);
    acc.ax(j) -=  dx*rr3*partc.mass(i);
    acc.ay(j) -=  dy*rr3*partc.mass(i);
    acc.az(j) -=  dz*rr3*partc.mass(i);
}
}*/
    template <typename Callable>
    class BasicEoM {
    public:
        template <typename Particles, typename Acc>
        void eval_acc(Particles& partc, Acc& acc){
            lambda_(partc,acc);
        }
    private:
        Callable lambda_;
    };

    class BasicEoM {
    public:
        template <typename Particles, typename Acc>
        void eval_acc(Particles& partc, Acc& acc){
            size_t num = partc.number();
            for(size_t i = 0 ; i < num; ++i){
                for(size_t j = i+1 ; j < num; ++j){
                    auto dx  = partc.px(j) - partc.px(i);
                    auto dy  = partc.py(j) - partc.py(i);
                    auto dz  = partc.pz(j) - partc.pz(i);
                    auto r   = sqrt(dx*dx+dy*dy+dz*dz);
                    auto rr3 = 1/(r*r*r);
                    acc.ax(i) +=  dx*rr3*partc.mass(j);
                    acc.ay(i) +=  dy*rr3*partc.mass(j);
                    acc.az(i) +=  dz*rr3*partc.mass(j);
                    acc.ax(j) -=  dx*rr3*partc.mass(i);
                    acc.ay(j) -=  dy*rr3*partc.mass(i);
                    acc.az(j) -=  dz*rr3*partc.mass(i);
                }
            }
        }
    };
}

#endif //SPACEHUB_EOM_H
