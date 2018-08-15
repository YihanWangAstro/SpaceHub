
#ifndef MACROS_H
#define MACROS_H

/*constexpr double PI        = 3.14159265358979323;
constexpr double AU        = (PI / 648000);
constexpr double PC        = 1;
constexpr double M_SOLAR   = 1;
constexpr double M_JUPITER = 0.9547919E-3;
constexpr double R_SOLAR   = 2.25461E-8;
constexpr double YEAR      = 6.694685210039141E-08;
constexpr double DAY       = YEAR / 365.25636042;
constexpr double G         = 1;
constexpr double V_UNIT    = 6.54589713446219E-2;
constexpr double C         = 299792.458 / V_UNIT;
constexpr double KM        = 3.2407557442395564e-14;*/

namespace SpaceH
{
    
    namespace Const
    {
        constexpr double PI        = 3.14159265358979323;
    }
    
    namespace Unit
    {
        constexpr double AU        = 1;
        constexpr double M_SOLAR   = 1;
        constexpr double YEAR      = 2*Const::PI;
        
        constexpr double GYR       = 1e9*YEAR;
        constexpr double MYR       = 1e6*YEAR;
        constexpr double MONTH     = YEAR / 12;
        constexpr double DAY       = YEAR / 365.25636042;
        constexpr double HUBBLETIME= 14.7*GYR;
        
        constexpr double KM        = 1 / 149597870.7;
        constexpr double PC        = 648000/Const::PI;
        constexpr double R_SOLAR   = 6.957E5*KM;
        
        constexpr double M_EARTH   = 3.003E-6*M_SOLAR;
        constexpr double M_MECURY  = 0.055*M_EARTH;
        constexpr double M_VENUS   = 0.815*M_EARTH;
        constexpr double M_MARS    = 0.107*M_EARTH;
        constexpr double M_JUPITER = 317.8*M_EARTH;
        constexpr double M_SATURN  = 95.16*M_EARTH;
        constexpr double M_NEPTUNE = 17.15*M_EARTH;
        constexpr double M_MOON    = 0.012300*M_EARTH;
        
        constexpr double V_UNIT    = 2.9784651272402163E1;
        constexpr double KMS       = 1 / V_UNIT;
    }
    
    namespace Const
    {
        constexpr double G         = 1;
        constexpr double C         = 299792.458 * Unit::KMS ;
    }
}
#endif
