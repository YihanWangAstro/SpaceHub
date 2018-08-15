#ifndef ORBITS_H
#define ORBITS_H
#include "../protoType.h"
#include "../ownMath.h"
#include "../macros.h"
#include <math.h>

namespace SpaceH
{
    namespace Orbits
    {
        template<typename Vector, typename Scalar>
        void eulerRotate(Vector& v, const Scalar phi, const Scalar theta, const Scalar psi)
        {
            Vector tmp(0.0, 0.0, 0.0);
            
            tmp.x = v.x*(cos(phi)*cos(psi) - sin(phi)*cos(theta)*sin(psi))
                  - v.y*(cos(phi)*sin(psi) + sin(phi)*cos(theta)*cos(psi))
                  + v.z*(sin(phi)*sin(theta));
            
            tmp.y = v.x*(sin(phi)*cos(psi) + cos(phi)*cos(theta)*sin(psi))
                  - v.y*(sin(phi)*sin(psi) - cos(phi)*cos(theta)*cos(psi))
                  - v.z*(cos(phi)*sin(theta));
            
            tmp.z = v.x*sin(theta)*sin(psi) + v.y*sin(theta)*cos(psi) + v.z*cos(theta);
            
            v = tmp;
        };
        
        template<typename Scalar>
        Scalar getRandomMeanAnomaly(Scalar e, Scalar Mmin, Scalar Mmax)
        {
            if(e >= 0)
                return SpaceH::uniform()*(Mmax-Mmin) + Mmin;
            else
            {
                SpaceH::errMsg("Eccentrcity cannot be negative, Nan or inf!", __FILE__, __LINE__);
                return 0;
            }
        }
        
        template<typename Scalar>
        Scalar getTrueAnomaly(Scalar E, Scalar e)
        {
            if(0 <= e && e < 1)
                return 2*atan2(sqrt(1+e)*sin(E*0.5), sqrt(1-e)*cos(0.5*E));
            else if(e > 1)
                return 2*atan2(sqrt(1+e)*sinh(E*0.5), sqrt(e-1)*cosh(0.5*E));
            else if(fabs(e - 1) < SpaceH::epsilon<Scalar>::value)
                return 2*atan(0.5*E);
            else
            {
                SpaceH::errMsg("Eccentrcity cannot be negative, Nan or inf!", __FILE__, __LINE__);
                return 0;
            }
        }
        
        template<typename Scalar>
        Scalar getEccentricAnomaly(Scalar M, Scalar e)
        {
            if(0 <= e && e < 1)//newton iteration may encounter stationary point
                return SpaceH::root_dichotom( [&](Scalar x)->Scalar {return x - e*sin(x) - M;} );//find this function in ownMath.h
            else if(e > 1)
                return SpaceH::root_dichotom( [&](Scalar x)->Scalar {return e*sinh(x) - x - M;} );
            else if(fabs(e - 1) < SpaceH::epsilon<Scalar>::value)
                return SpaceH::root_dichotom( [&](Scalar x)->Scalar {return x + x*x*x/3 - M;} );
            else
            {
                SpaceH::errMsg("Eccentrcity cannot be negative, Nan or inf!", __FILE__, __LINE__);
                return 0;
            }
        }
    }
    
    template<typename Dtype>
    struct Particle
    {
        /* Typedef */
        using type   = SpaceH::ProtoType<Dtype, SpaceH::DYNAMICAL>;
        using Scalar = typename type::Scalar;
        using Vector = typename type::Vector;
        /* Typedef */
        
        Vector pos;
        Vector vel;
        Scalar mass;
        Scalar radius;
        
        friend std::ostream& operator<<(std::ostream& os, const Particle& P)
        {
            os << P.pos/Unit::AU << ' ' << P.vel/Unit::KMS << ' ' << P.mass/Unit::M_SOLAR << ' ' << P.radius/Unit::AU;
            return os;
        }
    };
    
    template<typename Dtype>
    struct Kepler
    {
        /* Typedef */
        using type   = SpaceH::ProtoType<Dtype, SpaceH::DYNAMICAL>;
        using Scalar = typename type::Scalar;
        using Vector = typename type::Vector;
        /* Typedef */
        
        inline const Scalar a() const { return a_; }
        
        inline const Scalar b() const { return b_; }
        
        inline const Scalar p() const { return p_; }//semi-latus rectum
        
        inline const Scalar e() const { return e_; }
        
        inline const Scalar trueAnomaly() const { return trueAnomaly_; }
        
        inline const Scalar phi() const { return phi_; }
        
        inline const Scalar theta() const { return theta_; }
        
        inline const Scalar psi() const { return psi_; }
        
        Kepler() = delete;
        
        Kepler(Scalar m1, Scalar m2, Scalar p, Scalar e)
        {
            calcuOrbitalParameter(m1, m2, p, e);
            phi_         = SpaceH::uniform()*2*Const::PI - Const::PI;
            theta_       = acos(SpaceH::uniform()*2 - 1);
            psi_         = SpaceH::uniform()*2*Const::PI - Const::PI;
            Scalar M     = Orbits::getRandomMeanAnomaly(e, -Unit::HUBBLETIME/T_, Unit::HUBBLETIME/T_);
            Scalar E     = Orbits::getEccentricAnomaly(M, e_);
            trueAnomaly_ = Orbits::getTrueAnomaly(E, e_);
            createOrbit(m1, m2);
        }
        
        Kepler(Scalar m1, Scalar m2, Scalar p, Scalar e, Scalar phi, Scalar theta, Scalar psi, Scalar trueAnomaly = Const::PI)
        {
            calcuOrbitalParameter(m1, m2, p, e);
            phi_         = phi;
            theta_       = theta;
            psi_         = psi;
            trueAnomaly_ = trueAnomaly;
            createOrbit(m1, m2);
        }
        
        void moveOrbitTo(const Vector& newCMPos, const Vector& newCMVel)
        {
            moveToCentreMassCoord();
            P1.pos += newCMPos, P1.vel += newCMVel;
            P2.pos += newCMPos, P2.vel += newCMVel;
        }
        
        void moveToCentreMassCoord()
        {
            Vector CMP = getCentreMassPos();
            Vector CMV = getCentreMassVel();
            P1.pos -= CMP, P1.vel -= CMV;
            P2.pos -= CMP, P2.vel -= CMV;
        }
        
        void moveToPrimaryCentre()
        {
            if(P1.mass > P2.mass)
            {
                P2.pos -= P1.pos, P1.pos.setZero();
                P2.vel -= P1.vel, P1.vel.setZero();
            }
            else
            {
                P1.pos -= P2.pos, P2.pos.setZero();
                P1.vel -= P2.vel, P2.vel.setZero();
            }
        }
        
        inline Vector getCentreMassPos()
        {
            Scalar totalMass = P1.mass + P2.mass;
            return Vector((P1.mass*P1.pos + P2.mass*P2.pos)/totalMass);
        }
        
        inline Vector getCentreMassVel()
        {
            Scalar totalMass = P1.mass + P2.mass;
            return Vector((P1.mass*P1.vel + P2.mass*P2.vel)/totalMass);
        }
        
        inline const Particle<Dtype>& primary() const
        {
            if(P1.mass > P2.mass)
                return P1;
            else
                return P2;
        }
        
        inline const Particle<Dtype>& secondary() const
        {
            if(P1.mass < P2.mass)
                return P1;
            else
                return P2;
        }
        
        void randomPhi()
        {
            phi_   = SpaceH::uniform()*2*Const::PI - Const::PI;
            createOrbit(P1.mass, P2.mass);
        }
        
        void randomTheta()
        {
            theta_ = acos(SpaceH::uniform()*2 - 1);
            createOrbit(P1.mass, P2.mass);
        }
        
        void randomPsi()
        {
            psi_   = SpaceH::uniform()*2*Const::PI - Const::PI;
            createOrbit(P1.mass, P2.mass);
        }
        
        void randomAngles()
        {
            phi_   = SpaceH::uniform()*2*Const::PI - Const::PI;
            theta_ = acos(SpaceH::uniform()*2 - 1);
            psi_   = SpaceH::uniform()*2*Const::PI - Const::PI;
            createOrbit(P1.mass, P2.mass);
        }
        
        void randomPhase(Scalar pastTime = -Unit::HUBBLETIME, Scalar futureTime = Unit::HUBBLETIME)
        {
            Scalar M = Orbits::getRandomMeanAnomaly(e_, pastTime/T_, futureTime/T_);
            Scalar E = Orbits::getEccentricAnomaly(M, e_);
            trueAnomaly_ = Orbits::getTrueAnomaly(E, e_);
            createOrbit(P1.mass, P2.mass);
        }
    private:
        void calcuOrbitalParameter(Scalar m1, Scalar m2, Scalar p, Scalar e)
        {
            checkParameter(p, e);
            P1.mass = m1*Unit::M_SOLAR;
            P2.mass = m2*Unit::M_SOLAR;
            p_ = p*Unit::AU;
            e_ = e;
            if(fabs(e_ - 1) > SpaceH::epsilon<Scalar>::value)
            {
                a_ = p_/(1-e_*e_);
                b_ = fabs(a_)*sqrt(fabs(1-e_*e_));
                T_ = sqrt(fabs(a_*a_*a_/Const::G/(P1.mass+P2.mass)));
            }
            else
            {
                a_ = 1.0/0.0;
                b_ = 0.5*p_;
                T_ =sqrt(fabs(0.25*p_*p_*p_/Const::G/(P1.mass+P2.mass)));
            }
        }
        
        void checkParameter(Scalar p, Scalar e)
        {
            if(p < 0)
                SpaceH::errMsg("semi-latus rectum cannot be negative", __FILE__, __LINE__);
            if(e < 0)
                SpaceH::errMsg("Eccentrcity cannot be negative!", __FILE__, __LINE__);
        }
        
        void createOrbit(Scalar m1, Scalar m2)
        {
            Scalar u = (m1+m2)*Const::G;
            
            Scalar r = p_/(1 + e_*cos(trueAnomaly_)) ;
            P2.pos.x = r*cos(trueAnomaly_);
            P2.pos.y = r*sin(trueAnomaly_);
            P2.pos.z = 0;
            
            Scalar v  = sqrt(u/p_);
            P2.vel.x  = -v*sin(trueAnomaly_);
            P2.vel.y  = v*(e_ + cos(trueAnomaly_));
            P2.vel.z  = 0;
        
            P1.pos.setZero(), P1.vel.setZero();
            
            P1.mass = m1*Unit::M_SOLAR, P1.radius = 0;
            P2.mass = m2*Unit::M_SOLAR, P2.radius = 0;
            
            Orbits::eulerRotate(P2.pos, phi_, theta_, psi_);
            Orbits::eulerRotate(P2.vel, phi_, theta_, psi_);
            moveToCentreMassCoord();
        }
    private:
        Particle<Dtype> P1;
        Particle<Dtype> P2;
        
        Scalar e_;
        Scalar a_;
        Scalar b_;
        Scalar p_;
        Scalar T_;
        Scalar trueAnomaly_;
        Scalar phi_;
        Scalar theta_;
        Scalar psi_;
    };
}
#endif
