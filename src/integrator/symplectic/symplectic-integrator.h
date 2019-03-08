//
// Created by yihan on 3/8/19.
//

#ifndef SPACEHUB_SYMPLECTIC_INTEGRATOR_H
#define SPACEHUB_SYMPLECTIC_INTEGRATOR_H

#include "../../dev-tools.h"
#include "../../particle-system.h"

namespace SpaceH::integrator{

    template<typename Derived, size_t Order>
    class SymIntegrator{
    public:
        static constexpr int order{Order};

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
