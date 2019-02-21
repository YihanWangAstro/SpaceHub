
#ifndef SYMPLECTIC_6TH_INTEGRATOR_H
#define SYMPLECTIC_6TH_INTEGRATOR_H
#include "../../dev_tools.h"
/** @brief Sixth order symplectic integrator */
template<typename ParticSys>
class symplectic6th {
public:
    /* Typedef */
    SPACEHUB_USING_TYPE_SYSTEM_OF(ParticSys);
    /* Typedef */

    /*Template parameter check*/
    /*Template parameter check*/

    /** @brief Order of the integrator*/
    static constexpr int order{6};

    void integrate(ParticSys &particles, Scalar stepLength);
};

/** @brief Interface to integrate particle system
 *
 *  This function integrate the particle system for one step with DKD leapfrog second order symplectic algorithm.
 *  @param particles  Particle system need to be integrated.
 *  @param stepLength Step size for integration.
 */
template<typename ParticSys>
void symplectic6th<ParticSys>::integrate(ParticSys &particles, Scalar stepLength) {
    /*unroll loop manually*/
    particles.drift(3.9225680523877998E-1 * stepLength);
    particles.kick(7.8451361047755996E-1 * stepLength);
    particles.drift(5.1004341191845848E-1 * stepLength);
    particles.kick(2.3557321335935699E-1 * stepLength);
    particles.drift(-4.7105338540975655E-1 * stepLength);
    particles.kick(-1.1776799841788701E0 * stepLength);
    particles.drift(6.8753168252518093E-2 * stepLength);
    particles.kick(1.3151863206839063E0 * stepLength);
    particles.drift(6.8753168252518093E-2 * stepLength);
    particles.kick(-1.1776799841788701E0 * stepLength);
    particles.drift(-4.7105338540975655E-1 * stepLength);
    particles.kick(2.3557321335935699E-1 * stepLength);
    particles.drift(5.1004341191845848E-1 * stepLength);
    particles.kick(7.8451361047755996E-1 * stepLength);
    particles.drift(3.9225680523877998E-1 * stepLength);
}

#endif
