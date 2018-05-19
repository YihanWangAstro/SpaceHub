////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Filename:size_tegrator.h                                                                                               //
//Author:Yihan Wang                                                                                                   //
//                                                                                                                    //
//                                                                                                                    //
//Description:                                                                                                        //
//                                                                                                                    //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef SYMPLECTIC_2TH_INTEGRATOR_H
#define SYMPLECTIC_2TH_INTEGRATOR_H
/** @brief Second order symplectic integrator */
template <typename ParticSys>
class symplectic2th
{
    typedef typename ParticSys::Scalar Scalar;
public:
    /** @brief Order of the integrator*/
    static const int order{2};
    void integrate(ParticSys& particles, Scalar stepLength);
};
/** @brief Interface to integrate particle system
 *
 *  This function integrate the particle system for one step with DKD leapfrog second order symplectic algorithm.
 *  @param particles  Particle system need to be integrated.
 *  @param stepLength Step size for integration.
 */
template <typename ParticSys>
void symplectic2th<ParticSys>::integrate(ParticSys& particles, Scalar stepLength)
{
    particles.advancePos(0.5*stepLength);
    particles.advanceVel(stepLength);
    particles.advancePos(0.5*stepLength);
}
#endif
