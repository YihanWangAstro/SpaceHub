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
        auto iterate(T& particles, typename T::Scalar macro_step_size) -> typename T::Scalar {
            static_assert(is_particle_system<T>::value, "Passing non paritcle-system-type!");
            return static_cast<Derived*>(this)->impl_iterate(particles, macro_step_size);
        }
    private:
        OdeIterator() = default;
        friend Derived;
    };

    template <typename>
    struct is_ode_iterator : public std::false_type { };

    template <typename T>
    struct is_ode_iterator<OdeIterator<T>> : public std::true_type { };
}
#endif //SPACEHUB_ODE_ITERATOR_H
