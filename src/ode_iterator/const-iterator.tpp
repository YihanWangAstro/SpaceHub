
#ifndef CONSTITERATOR_H
#define CONSTITERATOR_H

#include "../dev-tools.h"
#include "ode-iterator.h"
namespace SpaceH::odeIterator{

    template<typename Integrator>
    class ConstOdeIterator : public OdeIterator<ConstOdeIterator<Integrator>> {
    public:
        template <typename T>
        auto impl_iterate(ParticleSystem<T>& particles, typename T::Scalar macroStepSize) -> typename T::Scalar {
            integrator_.integrate(particles, macroStepSize);
            return macroStepSize;
        }
    private:
        Integrator integrator_;
    };
}
#endif
