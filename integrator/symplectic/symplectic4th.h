////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Filename:integrator.h                                                                                               //
//Author:Yihan Wang                                                                                                   //
//                                                                                                                    //
//                                                                                                                    //
//Description:                                                                                                        //
//                                                                                                                    //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef SYMPLECTIC_4TH_INTEGRATOR_H
#define SYMPLECTIC_4TH_INTEGRATOR_H
/** @brief Fourth order symplectic integrator */
template <typename ParticSys>
class symplectic4th
{
public:
    /* Typedef */
    using type   = typename ParticSys::type;
    using Scalar = typename type::Scalar;
    /* Typedef */
    
    /*Template parameter check*/
    /*Template parameter check*/
    
    /** @brief Order of the integrator*/
    static const int order{4};
    void integrate(ParticSys& particles, Scalar stepLength);
};

/** @brief Interface to integrate particle system
 *
 *  This function integrate the particle system for one step with DKD leapfrog second order symplectic algorithm.
 *  @param particles  Particle system need to be integrated.
 *  @param stepLength Step size for integration.
 */
template <typename ParticSys>
void symplectic4th<ParticSys>::integrate(ParticSys& particles, Scalar stepLength)
{   /*unroll loop manually*/
    particles.drift(6.7560359597983000E-1 * stepLength);
    particles.kick(1.3512071919596600E0 * stepLength);
    particles.drift(-1.7560359597983000E-1 * stepLength);
    particles.kick(-1.7024143839193200E0 * stepLength);
    particles.drift(-1.7560359597983000E-1 * stepLength);
    particles.kick(1.3512071919596600E0 * stepLength);
    particles.drift(6.7560359597983000E-1 * stepLength);
}
#endif
