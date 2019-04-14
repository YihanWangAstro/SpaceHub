//
// Created by yihan on 3/8/19.
//

#ifndef SPACEHUB_SYMPLECTIC_INTEGRATOR_H
#define SPACEHUB_SYMPLECTIC_INTEGRATOR_H

#include "../../particle-system.h"

namespace space::integrator{

    template<typename Derived>
    class SymIntegrator{
    public:
        static constexpr size_t order{Derived::order};

        template <typename T>
        void integrate(T& ptc, typename T::Scalar stepSize) {
            static_assert(is_particle_system_v<T>, "Passing non paritcle-system-type!");
            static_cast<Derived*>(this)->impl_integrate(ptc, stepSize);
        }

    private:
        SymIntegrator() = default;
        friend Derived;
    };

    template <typename T>
    constexpr bool is_sym_integrator_v = std::is_base_of_v<SymIntegrator<T>, T>;
}
#endif //SPACEHUB_SYMPLECTIC_INTEGRATOR_H
