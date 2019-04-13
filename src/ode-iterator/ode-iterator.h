//
// Created by yihan on 3/8/19.
//

#ifndef SPACEHUB_ODE_ITERATOR_H
#define SPACEHUB_ODE_ITERATOR_H

#include "../particle-system.h"

namespace space::odeIterator{
    template <typename Derived>
    class OdeIterator{
    public:
        template <typename T>
        auto iterate(T& particles, typename T::Scalar macro_step_size) -> typename T::Scalar {
            static_assert(is_particle_system_v<T>, "Passing non paritcle-system-type!");
            return static_cast<Derived*>(this)->impl_iterate(particles, macro_step_size);
        }
    private:
        OdeIterator() = default;
        friend Derived;
    };

    template <typename T>
    constexpr bool is_ode_iterator_v = std::is_base_of_v<OdeIterator<T>, T>;
}
#endif //SPACEHUB_ODE_ITERATOR_H
