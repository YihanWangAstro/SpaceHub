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
        void integrate(ParticleSystem<T>& ptc, typename T::Scalar stepSize) {
            static_cast<Derived*>(this)->impl_integrate(ptc, stepSize);
        }

    private:
        SymIntegrator() = default;
        friend Derived;
    };
}
#endif //SPACEHUB_SYMPLECTIC_INTEGRATOR_H
