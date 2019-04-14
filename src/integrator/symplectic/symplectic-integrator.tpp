
#ifndef SYMPLECTIC_2TH_INTEGRATOR_H
#define SYMPLECTIC_2TH_INTEGRATOR_H

#include "symplectic-integrator.h"

namespace space::integrator {

    class symplectic2nd : public SymIntegrator<symplectic2nd> {
    public:
        static constexpr size_t order{2};
        
        template<typename T>
        void impl_integrate(ParticleSystem<T> &particles, typename T::Scalar step_size) {
            particles.drift(0.5 * step_size);
            particles.kick(step_size);
            particles.drift(0.5 * step_size);
        }
    };

    class symplectic4th : public SymIntegrator<symplectic4th> {
    public:
        static constexpr size_t order{4};

        template<typename T>
        void impl_integrate(ParticleSystem<T> &particles, typename T::Scalar step_size) {
            particles.drift(6.7560359597983000E-1 * step_size);
            particles.kick(1.3512071919596600E0 * step_size);
            particles.drift(-1.7560359597983000E-1 * step_size);
            particles.kick(-1.7024143839193200E0 * step_size);
            particles.drift(-1.7560359597983000E-1 * step_size);
            particles.kick(1.3512071919596600E0 * step_size);
            particles.drift(6.7560359597983000E-1 * step_size);
        }
    };

    class symplectic6th : public SymIntegrator<symplectic6th> {
    public:
        static constexpr size_t order{6};

        template<typename T>
        void impl_integrate(ParticleSystem<T> &particles, typename T::Scalar step_size) {
            /*unroll loop manually*/
            particles.drift(3.9225680523877998E-1 * step_size);
            particles.kick(7.8451361047755996E-1 * step_size);
            particles.drift(5.1004341191845848E-1 * step_size);
            particles.kick(2.3557321335935699E-1 * step_size);
            particles.drift(-4.7105338540975655E-1 * step_size);
            particles.kick(-1.1776799841788701E0 * step_size);
            particles.drift(6.8753168252518093E-2 * step_size);
            particles.kick(1.3151863206839063E0 * step_size);
            particles.drift(6.8753168252518093E-2 * step_size);
            particles.kick(-1.1776799841788701E0 * step_size);
            particles.drift(-4.7105338540975655E-1 * step_size);
            particles.kick(2.3557321335935699E-1 * step_size);
            particles.drift(5.1004341191845848E-1 * step_size);
            particles.kick(7.8451361047755996E-1 * step_size);
            particles.drift(3.9225680523877998E-1 * step_size);
        }
    };

    class symplectic8th : public SymIntegrator<symplectic8th> {
    public:
        static constexpr size_t order{8};

        template<typename T>
        void impl_integrate(ParticleSystem<T> &particles, typename T::Scalar step_size) {
            /*unroll loop manually*/
            particles.drift(5.21213104349955048E-1 * step_size);
            particles.kick(1.04242620869991010E0 * step_size);
            particles.drift(1.43131625920352512E0 * step_size);
            particles.kick(1.82020630970713992E0 * step_size);
            particles.drift(9.88973118915378424E-1 * step_size);
            particles.kick(1.57739928123617007E-1 * step_size);
            particles.drift(1.29888362714548355E0 * step_size);
            particles.kick(2.44002732616735019E0 * step_size);
            particles.drift(1.21642871598513458E0 * step_size);
            particles.kick(-7.16989419708119989E-3 * step_size);
            particles.drift(-1.22708085895116059E0 * step_size);
            particles.kick(-2.44699182370524015E0 * step_size);
            particles.drift(-2.03140778260310517E0 * step_size);
            particles.kick(-1.61582374150096997E0 * step_size);
            particles.drift(-1.69832618404521085E0 * step_size);
            particles.kick(-1.78082862658945151E0 * step_size);
            particles.drift(-1.69832618404521085E0 * step_size);
            particles.kick(-1.61582374150096997E0 * step_size);
            particles.drift(-2.03140778260310517E0 * step_size);
            particles.kick(-2.44699182370524015E0 * step_size);
            particles.drift(-1.22708085895116059E0 * step_size);
            particles.kick(-7.16989419708119989E-3 * step_size);
            particles.drift(1.21642871598513458E0 * step_size);
            particles.kick(2.44002732616735019E0 * step_size);
            particles.drift(1.29888362714548355E0 * step_size);
            particles.kick(1.57739928123617007E-1 * step_size);
            particles.drift(9.88973118915378424E-1 * step_size);
            particles.kick(1.82020630970713992E0 * step_size);
            particles.drift(1.43131625920352512E0 * step_size);
            particles.kick(1.04242620869991010E0 * step_size);
            particles.drift(5.21213104349955048E-1 * step_size);
        }
    };

    class symplectic10th : public SymIntegrator<symplectic10th> {
    public:
        static constexpr size_t order{10};

        template<typename T>
        void impl_integrate(ParticleSystem<T> &particles, typename T::Scalar step_size) {
            /*unroll loop manually*/
            particles.drift(3.0610967201933609e-01 * step_size);
            particles.kick(6.1221934403867218e-01 * step_size);
            particles.drift(-9.4012698954724694e-02 * step_size);
            particles.kick(-8.0024474194812156e-01 * step_size);
            particles.drift(-6.6002635995076209e-01 * step_size);
            particles.kick(-5.1980797795340250e-01 * step_size);
            particles.drift(-1.5240397828727220e-01 * step_size);
            particles.kick(2.1500002137885812e-01 * step_size);
            particles.drift(-1.1750569210727700e-01 * step_size);
            particles.kick(-4.5001140559341213e-01 * step_size);
            particles.drift(2.2250778443570857e-01 * step_size);
            particles.kick(8.9502697446482926e-01 * step_size);
            particles.drift(5.1288848042847668e-01 * step_size);
            particles.kick(1.3074998639212410e-01 * step_size);
            particles.drift(3.3095796002497074e-01 * step_size);
            particles.kick(5.3116593365781739e-01 * step_size);
            particles.drift(-6.0050191119721985e-02 * step_size);
            particles.kick(-6.5126631589726136e-01 * step_size);
            particles.drift(-7.6956706144236287e-01 * step_size);
            particles.kick(-8.8786780698746448e-01 * step_size);
            particles.drift(-7.6872229417056015e-02 * step_size);
            particles.kick(7.3412334815335245e-01 * step_size);
            particles.drift(4.2477286784491525e-01 * step_size);
            particles.kick(1.1542238753647800e-01 * step_size);
            particles.drift(4.3160892192959932e-01 * step_size);
            particles.kick(7.4779545632272060e-01 * step_size);
            particles.drift(5.5434862753225678e-02 * step_size);
            particles.kick(-6.3692573081626924e-01 * step_size);
            particles.drift(-1.9288621063874828e-01 * step_size);
            particles.kick(2.5115330953877268e-01 * step_size);
            particles.drift(3.3904387248169282e-01 * step_size);
            particles.kick(4.2693443542461296e-01 * step_size);
            particles.drift(3.3904387248169282e-01 * step_size);
            particles.kick(2.5115330953877268e-01 * step_size);
            particles.drift(-1.9288621063874828e-01 * step_size);
            particles.kick(-6.3692573081626924e-01 * step_size);
            particles.drift(5.5434862753225678e-02 * step_size);
            particles.kick(7.4779545632272060e-01 * step_size);
            particles.drift(4.3160892192959932e-01 * step_size);
            particles.kick(1.1542238753647800e-01 * step_size);
            particles.drift(4.2477286784491525e-01 * step_size);
            particles.kick(7.3412334815335245e-01 * step_size);
            particles.drift(-7.6872229417056015e-02 * step_size);
            particles.kick(-8.8786780698746448e-01 * step_size);
            particles.drift(-7.6956706144236287e-01 * step_size);
            particles.kick(-6.5126631589726136e-01 * step_size);
            particles.drift(-6.0050191119721985e-02 * step_size);
            particles.kick(5.3116593365781739e-01 * step_size);
            particles.drift(3.3095796002497074e-01 * step_size);
            particles.kick(1.3074998639212410e-01 * step_size);
            particles.drift(5.1288848042847668e-01 * step_size);
            particles.kick(8.9502697446482926e-01 * step_size);
            particles.drift(2.2250778443570857e-01 * step_size);
            particles.kick(-4.5001140559341213e-01 * step_size);
            particles.drift(-1.1750569210727700e-01 * step_size);
            particles.kick(2.1500002137885812e-01 * step_size);
            particles.drift(-1.5240397828727220e-01 * step_size);
            particles.kick(-5.1980797795340250e-01 * step_size);
            particles.drift(-6.6002635995076209e-01 * step_size);
            particles.kick(-8.0024474194812156e-01 * step_size);
            particles.drift(-9.4012698954724694e-02 * step_size);
            particles.kick(6.1221934403867218e-01 * step_size);
            particles.drift(3.0610967201933609e-01 * step_size);
        }
    };
}
#endif
