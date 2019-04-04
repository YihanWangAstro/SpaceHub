//
// Created by yihan on 3/8/19.
//

#ifndef SPACEHUB_ODE_ITERATOR_H
#define SPACEHUB_ODE_ITERATOR_H

#include "../particle-system.h"

namespace SpaceH::OdeIterator{
    template <typename Derived>
    class OdeIterator{
    public:
        template <typename T>
        auto iterate(ParticleSystem<T>& particles, typename T::Scalar macro_step_size) -> typename T::Scalar {
            return static_cast<Derived*>(this)->impl_iterate(particles, macro_step_size);
        }
    private:
        OdeIterator() = default;
        friend Derived;
    };
}
#endif //SPACEHUB_ODE_ITERATOR_H
