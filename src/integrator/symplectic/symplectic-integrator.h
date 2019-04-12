//
// Created by yihan on 3/8/19.
//

#ifndef SPACEHUB_SYMPLECTIC_INTEGRATOR_H
#define SPACEHUB_SYMPLECTIC_INTEGRATOR_H

#include "../../particle-system.h"

namespace SpaceH::Integrator{

    template<typename Derived, size_t Order>
    class SymIntegrator{
    public:
        static constexpr size_t order{Order};

        template <typename T>
        void integrate(T& ptc, typename T::Scalar stepSize) {
            static_assert(is_particle_system<T>::value, "Passing non paritcle-system-type!");
            static_cast<Derived*>(this)->impl_integrate(ptc, stepSize);
        }

    private:
        SymIntegrator() = default;
        friend Derived;
    };

    template <typename>
    struct is_sym_integrator : public std::false_type { };

    template <typename T, size_t Order>
    struct is_sym_integrator<SymIntegrator<T, Order>> : public std::true_type { };
}
#endif //SPACEHUB_SYMPLECTIC_INTEGRATOR_H
