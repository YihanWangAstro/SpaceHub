
#ifndef CONSTITERATOR_H
#define CONSTITERATOR_H

#include "../dev-tools.hpp"
#include "ode-iterator.hpp"
namespace space::odeIterator{

    template<typename Integrator>
    class ConstOdeIterator : public OdeIterator<ConstOdeIterator<Integrator>> {
    public:
        template <typename T>
        auto impl_iterate(T& particles, typename T::Scalar macro_step_size) -> typename T::Scalar {
            static_assert(is_particle_system_v<T>, "Passing non paritcle-system-type!");
            integrator_.integrate(particles, macro_step_size);
            return macro_step_size;
        }
    private:
        Integrator integrator_;
    };
}
#endif
