////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Filename:integrator.h                                                                                               //
//Author:Yihan Wang                                                                                                   //
//                                                                                                                    //
//                                                                                                                    //
//Description:                                                                                                        //
//                                                                                                                    //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef SYMPLECTIC_10TH_INTEGRATOR_H
#define SYMPLECTIC_10TH_INTEGRATOR_H
/** @brief Eighth order symplectic integrator */
template <typename ParticSys>
class symplectic10th
{
public:
    using type = typename ParticSys::type;
    using Scalar = typename type::Scalar;
    
    /** @brief Order of the integrator*/
    static const int order{10};
    void integrate(ParticSys& particles, Scalar stepLength);
};

/** @brief Interface to integrate particle system
 *
 *  This function integrate the particle system for one step with DKD leapfrog second order symplectic algorithm.
 *  @param particles  Particle system need to be integrated.
 *  @param stepLength Step size for integration.
 */
template <typename ParticSys>
void symplectic10th<ParticSys>::integrate(ParticSys& particles, Scalar stepLength)
{
    particles.drift(3.0610967201933609e-01 * stepLength);
    particles.kick(6.1221934403867218e-01 * stepLength);
    particles.drift(-9.4012698954724694e-02 * stepLength);
    particles.kick(-8.0024474194812156e-01 * stepLength);
    particles.drift(-6.6002635995076209e-01 * stepLength);
    particles.kick(-5.1980797795340250e-01 * stepLength);
    particles.drift(-1.5240397828727220e-01 * stepLength);
    particles.kick(2.1500002137885812e-01 * stepLength);
    particles.drift(-1.1750569210727700e-01 * stepLength);
    particles.kick(-4.5001140559341213e-01 * stepLength);
    particles.drift(2.2250778443570857e-01 * stepLength);
    particles.kick(8.9502697446482926e-01 * stepLength);
    particles.drift(5.1288848042847668e-01 * stepLength);
    particles.kick(1.3074998639212410e-01 * stepLength);
    particles.drift(3.3095796002497074e-01 * stepLength);
    particles.kick(5.3116593365781739e-01 * stepLength);
    particles.drift(-6.0050191119721985e-02 * stepLength);
    particles.kick(-6.5126631589726136e-01 * stepLength);
    particles.drift(-7.6956706144236287e-01 * stepLength);
    particles.kick(-8.8786780698746448e-01 * stepLength);
    particles.drift(-7.6872229417056015e-02 * stepLength);
    particles.kick(7.3412334815335245e-01 * stepLength);
    particles.drift(4.2477286784491525e-01 * stepLength);
    particles.kick(1.1542238753647800e-01 * stepLength);
    particles.drift(4.3160892192959932e-01 * stepLength);
    particles.kick(7.4779545632272060e-01 * stepLength);
    particles.drift(5.5434862753225678e-02 * stepLength);
    particles.kick(-6.3692573081626924e-01* stepLength);
    particles.drift(-1.9288621063874828e-01 * stepLength);
    particles.kick(2.5115330953877268e-01 * stepLength);
    particles.drift(3.3904387248169282e-01 * stepLength);
    particles.kick(4.2693443542461296e-01 * stepLength);
    particles.drift(3.3904387248169282e-01 * stepLength);
    particles.kick(2.5115330953877268e-01 * stepLength);
    particles.drift(-1.9288621063874828e-01 * stepLength);
    particles.kick(-6.3692573081626924e-01 * stepLength);
    particles.drift(5.5434862753225678e-02 * stepLength);
    particles.kick(7.4779545632272060e-01 * stepLength);
    particles.drift(4.3160892192959932e-01 * stepLength);
    particles.kick(1.1542238753647800e-01 * stepLength);
    particles.drift(4.2477286784491525e-01 * stepLength);
    particles.kick(7.3412334815335245e-01 * stepLength);
    particles.drift(-7.6872229417056015e-02 * stepLength);
    particles.kick(-8.8786780698746448e-01 * stepLength);
    particles.drift(-7.6956706144236287e-01 * stepLength);
    particles.kick(-6.5126631589726136e-01 * stepLength);
    particles.drift(-6.0050191119721985e-02 * stepLength);
    particles.kick(5.3116593365781739e-01 * stepLength);
    particles.drift(3.3095796002497074e-01 * stepLength);
    particles.kick(1.3074998639212410e-01 * stepLength);
    particles.drift(5.1288848042847668e-01 * stepLength);
    particles.kick(8.9502697446482926e-01 * stepLength);
    particles.drift(2.2250778443570857e-01 * stepLength);
    particles.kick(-4.5001140559341213e-01 * stepLength);
    particles.drift(-1.1750569210727700e-01 * stepLength);
    particles.kick(2.1500002137885812e-01 * stepLength);
    particles.drift(-1.5240397828727220e-01 * stepLength);
    particles.kick(-5.1980797795340250e-01 * stepLength);
    particles.drift(-6.6002635995076209e-01 * stepLength);
    particles.kick(-8.0024474194812156e-01 * stepLength);
    particles.drift(-9.4012698954724694e-02 * stepLength);
    particles.kick(6.1221934403867218e-01 * stepLength);
    particles.drift(3.0610967201933609e-01 * stepLength);
}
#endif
