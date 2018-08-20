#include <fstream>
#include "../../initCreator/orbits.h"
#include "../../initCreator/InitFileCreator.h"

using namespace SpaceH;
using Scalar    = double;
using particle  = Particle<Scalar>;
using Creator   = InitCreator<particle>;
using KObt      = Kepler<Scalar>;

constexpr size_t N = 5000;
constexpr size_t orbitsNum = 10;
constexpr size_t interval = N/orbitsNum;
int main()
{
    {
        Scalar e = 0.0167086;
        Scalar p = fabs(1-e*e)*1.000001018* Unit::AU;
        
        KObt orbit(Unit::M_SOLAR, Unit::M_EARTH, p ,e, 0, 0, 0);
        Creator crt;
        
        crt.addParticle(orbit.primary());
        crt.addParticle(orbit.secondary());
        crt.writeToFile("solar_earth.init");
    }
    
    {
        Scalar e = 0;
        Scalar p = fabs(1-e*e)* Unit::AU;
        
        KObt orbit(Unit::M_SOLAR, Unit::M_EARTH, p ,e, 0, 0, 0);
        Creator crt;
        
        crt.addParticle(orbit.primary());
        crt.addParticle(orbit.secondary());
        crt.writeToFile("circular.init");
    }
    
    {
        Scalar e = 0.999;
        Scalar p = fabs(1-e*e)* Unit::AU;
        
        KObt orbit(Unit::M_SOLAR, Unit::M_EARTH, p ,e, 0, 0, 0);
        Creator crt;
        
        crt.addParticle(orbit.primary());
        crt.addParticle(orbit.secondary());
        crt.writeToFile("elliptic.init");
    }
    
    {
        Scalar e = 0.5;
        Scalar p = fabs(1-e*e)* Unit::AU;
        KObt orbit1(1*Unit::M_SOLAR, Unit::M_EARTH, p, e, 0, 0, 0);
        
        e = 0.01;
        p = 0.001 * fabs(1-e*e) * Unit::AU;
        KObt orbit2(Unit::M_EARTH, Unit::M_EARTH, p, e, 0, 0.5*Const::PI, 0);

        orbit2.moveOrbitTo(orbit1.secondary());
        
        Creator crt;
        
        crt.addParticle(orbit1.primary());
        crt.addParticle(orbit2.primary());
        crt.addParticle(orbit2.secondary());
        crt.writeToFile("Kozai.init");
    }
    
    
}
