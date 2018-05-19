////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Filename:integrator.h                                                                                               //
//Author:Yihan Wang                                                                                                   //
//                                                                                                                    //
//                                                                                                                    //
//Description:                                                                                                        //
//                                                                                                                    //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef SYMPLECTIC_6TH_INTEGRATOR_H
#define SYMPLECTIC_6TH_INTEGRATOR_H
/** @brief Sixth order symplectic integrator */
template <typename ParticSys>
class symplectic6th
{
    //////////////////////////////////Interface/////////////////////////////////////
public:
    /** @brief Order of the integrator*/
   static const int order{6};
private:
    void integrate(ParticSys& particles, double stepLength);
};
/** @brief Interface to integrate particle system
 *
 *  This function integrate the particle system for one step with DKD leapfrog second order symplectic algorithm.
 *  @param particles  Particle system need to be integrated.
 *  @param stepLength Step size for integration.
 */
template <typename ParticSys>
void symplectic6th<ParticSys>::integrate(ParticSys& particles, double stepLength)
{
    particles.AdvancePos(3.9225680523877998E-1*stepLength);
    particles.AdvanceVel(7.8451361047755996E-1*stepLength);
    particles.AdvancePos(5.1004341191845848E-1*stepLength);
    particles.AdvanceVel(2.3557321335935699E-1*stepLength);
    particles.AdvancePos(-4.7105338540975655E-1*stepLength);
    particles.AdvanceVel(-1.1776799841788701E0*stepLength);
    particles.AdvancePos(6.8753168252518093E-2*stepLength);
    particles.AdvanceVel(1.3151863206839063E0*stepLength);
    particles.AdvancePos(6.8753168252518093E-2*stepLength);
    particles.AdvanceVel(-1.1776799841788701E0*stepLength);
    particles.AdvancePos(-4.7105338540975655E-1*stepLength);
    particles.AdvanceVel(2.3557321335935699E-1*stepLength);
    particles.AdvancePos(5.1004341191845848E-1*stepLength);
    particles.AdvanceVel(7.8451361047755996E-1*stepLength);
    particles.AdvancePos(3.9225680523877998E-1*stepLength);
}
#endif
