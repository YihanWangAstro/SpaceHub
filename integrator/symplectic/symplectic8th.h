////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Filename:integrator.h                                                                                               //
//Author:Yihan Wang                                                                                                   //
//                                                                                                                    //
//                                                                                                                    //
//Description:                                                                                                        //
//                                                                                                                    //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef SYMPLECTIC_8TH_INTEGRATOR_H
#define SYMPLECTIC_8TH_INTEGRATOR_H
/** @brief Eighth order symplectic integrator */
template <typename ParticSys>
class symplectic8th
{
    
public:
    /** @brief Order of the integrator*/
    static const int order{8};
    void integrate(ParticSys& particles, double stepLength);
};

/** @brief Interface to integrate particle system
 *
 *  This function integrate the particle system for one step with DKD leapfrog second order symplectic algorithm.
 *  @param particles  Particle system need to be integrated.
 *  @param stepLength Step size for integration.
 */
template <typename ParticSys>
void symplectic8th<ParticSys>::integrate(ParticSys& particles, double stepLength)
{   /*unroll loop*/
    particles.advancePos(5.21213104349955048E-1 * stepLength);
    particles.advanceVel(1.04242620869991010E0 * stepLength);
    particles.advancePos(1.43131625920352512E0 * stepLength);
    particles.advanceVel(1.82020630970713992E0 * stepLength);
    particles.advancePos(9.88973118915378424E-1 * stepLength);
    particles.advanceVel(1.57739928123617007E-1 * stepLength);
    particles.advancePos(1.29888362714548355E0 * stepLength);
    particles.advanceVel(2.44002732616735019E0 * stepLength);
    particles.advancePos(1.21642871598513458E0 * stepLength);
    particles.advanceVel(-7.16989419708119989E-3 * stepLength);
    particles.advancePos(-1.22708085895116059E0 * stepLength);
    particles.advanceVel(-2.44699182370524015E0 * stepLength);
    particles.advancePos(-2.03140778260310517E0 * stepLength);
    particles.advanceVel(-1.61582374150096997E0 * stepLength);
    particles.advancePos(-1.69832618404521085E0 * stepLength);
    particles.advanceVel(-1.78082862658945151E0 * stepLength);
    particles.advancePos(-1.69832618404521085E0 * stepLength);
    particles.advanceVel(-1.61582374150096997E0 * stepLength);
    particles.advancePos(-2.03140778260310517E0 * stepLength);
    particles.advanceVel(-2.44699182370524015E0 * stepLength);
    particles.advancePos(-1.22708085895116059E0 * stepLength);
    particles.advanceVel(-7.16989419708119989E-3 * stepLength);
    particles.advancePos(1.21642871598513458E0 * stepLength);
    particles.advanceVel(2.44002732616735019E0 * stepLength);
    particles.advancePos(1.29888362714548355E0 * stepLength);
    particles.advanceVel(1.57739928123617007E-1 * stepLength);
    particles.advancePos(9.88973118915378424E-1 * stepLength);
    particles.advanceVel(1.82020630970713992E0 * stepLength);
    particles.advancePos(1.43131625920352512E0 * stepLength);
    particles.advanceVel(1.04242620869991010E0 * stepLength);
    particles.advancePos(5.21213104349955048E-1 * stepLength);
}
#endif
