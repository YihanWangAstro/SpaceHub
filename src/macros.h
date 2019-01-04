
#ifndef MACROS_H
#define MACROS_H

namespace SpaceH::Const
{
    constexpr double PI        = 3.14159265358979323;
}

namespace SpaceH::Unit
{
    constexpr double AU        = 1;
    constexpr double M_SOLAR   = 1;
    constexpr double YEAR      = 2*Const::PI;

    constexpr double KYR       = 1e3*YEAR;
    constexpr double MYR       = 1e6*YEAR;
    constexpr double GYR       = 1e9*YEAR;
    constexpr double MONTH     = YEAR / 12;
    constexpr double DAY       = YEAR / 365.25636042;
    constexpr double HOUR      = DAY / 24;
    constexpr double MINUTE    = HOUR / 60;
    constexpr double SECOND    = MINUTE /60;
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

    constexpr double DEG       = Const::PI/180;

    struct Fmt {
        double L;
        double M;
        double T;
        double V;
    };
    constexpr Fmt STD_UNIT = {AU, M_SOLAR, YEAR, KMS};
    constexpr Fmt UNITY_UNIT = {1.0, 1.0, 1.0, 1.0};
}

namespace SpaceH::Const
{
    constexpr double G         = 1;
    constexpr double C         = 299792.458 * Unit::KMS ;
}



#endif
