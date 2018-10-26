
#ifndef SYMPLECTIC_2TH_INTEGRATOR_H
#define SYMPLECTIC_2TH_INTEGRATOR_H
#include "../../dev_tools.h"
namespace SpaceH {
    /** @brief Second order symplectic integrator */
    template<typename ParticSys>
    class symplectic2th {
    public:
        /* Typedef */
        SPACEHUB_USING_TYPE_SYSTEM_OF(ParticSys);
        /* Typedef */

        /*Template parameter check*/
        /*Template parameter check*/

        /** @brief Order of the integrator*/
        static const int order{2};

        void integrate(ParticSys &particles, Scalar stepLength);
    };

    /** @brief Interface to integrate particle system
     *
     *  This function integrate the particle system for one step with DKD leapfrog second order symplectic algorithm.
     *  @param particles  Particle system need to be integrated.
     *  @param stepLength Step size for integration.
     */
    template<typename ParticSys>
    void symplectic2th<ParticSys>::integrate(ParticSys &particles, Scalar stepLength) {
        particles.drift(0.5 * stepLength);
        particles.kick(stepLength);
        particles.drift(0.5 * stepLength);
    }
}
#endif
