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
    using type = typename ParticSys::type;
    using Scalar = typename type::Scalar;
    
    /** @brief Order of the integrator*/
    static const int order{8};
    void integrate(ParticSys& particles, Scalar stepLength);
};

/** @brief Interface to integrate particle system
 *
 *  This function integrate the particle system for one step with DKD leapfrog second order symplectic algorithm.
 *  @param particles  Particle system need to be integrated.
 *  @param stepLength Step size for integration.
 */
template <typename ParticSys>
void symplectic8th<ParticSys>::integrate(ParticSys& particles, Scalar stepLength)
{   /*unroll loop*/
    particles.drift(5.21213104349955048E-1 * stepLength);
    particles.kick(1.04242620869991010E0 * stepLength);
    particles.drift(1.43131625920352512E0 * stepLength);
    particles.kick(1.82020630970713992E0 * stepLength);
    particles.drift(9.88973118915378424E-1 * stepLength);
    particles.kick(1.57739928123617007E-1 * stepLength);
    particles.drift(1.29888362714548355E0 * stepLength);
    particles.kick(2.44002732616735019E0 * stepLength);
    particles.drift(1.21642871598513458E0 * stepLength);
    particles.kick(-7.16989419708119989E-3 * stepLength);
    particles.drift(-1.22708085895116059E0 * stepLength);
    particles.kick(-2.44699182370524015E0 * stepLength);
    particles.drift(-2.03140778260310517E0 * stepLength);
    particles.kick(-1.61582374150096997E0 * stepLength);
    particles.drift(-1.69832618404521085E0 * stepLength);
    particles.kick(-1.78082862658945151E0 * stepLength);
    particles.drift(-1.69832618404521085E0 * stepLength);
    particles.kick(-1.61582374150096997E0 * stepLength);
    particles.drift(-2.03140778260310517E0 * stepLength);
    particles.kick(-2.44699182370524015E0 * stepLength);
    particles.drift(-1.22708085895116059E0 * stepLength);
    particles.kick(-7.16989419708119989E-3 * stepLength);
    particles.drift(1.21642871598513458E0 * stepLength);
    particles.kick(2.44002732616735019E0 * stepLength);
    particles.drift(1.29888362714548355E0 * stepLength);
    particles.kick(1.57739928123617007E-1 * stepLength);
    particles.drift(9.88973118915378424E-1 * stepLength);
    particles.kick(1.82020630970713992E0 * stepLength);
    particles.drift(1.43131625920352512E0 * stepLength);
    particles.kick(1.04242620869991010E0 * stepLength);
    particles.drift(5.21213104349955048E-1 * stepLength);
}
#endif
