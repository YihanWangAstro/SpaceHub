
#ifndef SPACEHUB_MACROS_H
#define SPACEHUB_MACROS_H

namespace space::consts {
    constexpr double pi = 3.14159265358979323;
}

namespace space::unit {
    constexpr double au = 1;
    constexpr double m_solar = 1;
    constexpr double year = 2 * consts::pi;

    constexpr double kyr = 1e3 * year;
    constexpr double Myr = 1e6 * year;
    constexpr double Gyr = 1e9 * year;
    constexpr double month = year / 12;
    constexpr double day = year / 365.25636042;
    constexpr double hr = day / 24;
    constexpr double min = hr / 60;
    constexpr double sec = min / 60;
    constexpr double hubble_t = 13.7 * Gyr;

    constexpr double km = 1 / 149597870.7;
    constexpr double pc = 648000 / consts::pi;
    constexpr double r_solar = 6.957E5 * km;

    constexpr double m_earth = 3.003E-6 * m_solar;
    constexpr double m_mercury = 0.055 * m_earth;
    constexpr double m_venus = 0.815 * m_earth;
    constexpr double m_mars = 0.107 * m_earth;
    constexpr double m_jupiter = 317.8 * m_earth;
    constexpr double m_saturn = 95.16 * m_earth;
    constexpr double m_neptune = 17.15 * m_earth;
    constexpr double m_moon = 0.012300 * m_earth;

    constexpr double v_unit = 2.9784651272402163E1;
    constexpr double kms = 1 / v_unit;

    constexpr double deg = consts::pi / 180.0;

}

namespace space::consts {
    constexpr double G = 1;
    constexpr double C = 299792.458 * unit::kms;
}


#endif
