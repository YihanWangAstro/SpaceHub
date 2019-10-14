
#ifndef SPACEHUB_MACROS_HPP
#define SPACEHUB_MACROS_HPP

/**
 * @namespace space::consts
 * Documentation for space
 */
namespace space::consts {
  constexpr double pi = 3.14159265358979323;
}

/**
 * @namespace space::unit
 * Documentation for space
 */
namespace space::unit {
  constexpr double au = 1;
  constexpr double m_solar = 1;
  constexpr double year = 2 * consts::pi;
  constexpr double deg = consts::pi / 180.0;

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
  constexpr double r_jupiter = 69911 * km;
  constexpr double r_neptune = 24622 * km;

  constexpr double v_unit = 2.9784651272402163E1;
  constexpr double kms = 1 / v_unit;

  constexpr double m_earth = 3.003E-6 * m_solar;
  constexpr double m_mercury = 0.055 * m_earth;
  constexpr double m_venus = 0.815 * m_earth;
  constexpr double m_mars = 0.107 * m_earth;
  constexpr double m_jupiter = 317.8 * m_earth;
  constexpr double m_saturn = 95.16 * m_earth;
  constexpr double m_neptune = 17.15 * m_earth;
  constexpr double m_moon = 0.012300 * m_earth;
}  // namespace space::unit

namespace space::consts {
  constexpr double a_jupiter = 5.2044 * unit::au;
  constexpr double e_jupiter = 0.0489;
  constexpr double LoAN_jupiter = 100.464 * unit::deg;
  constexpr double AoP_jupiter = 273.867 * unit::deg;
  constexpr double i_jupiter = 6.09 * unit::deg;

  constexpr double a_neptune = 30.11 * unit::au;
  constexpr double e_neptune = 0.009456;
  constexpr double LoAN_neptune = 131.784 * unit::deg;
  constexpr double AoP_neptune = 276.336 * unit::deg;
  constexpr double i_neptune = 6.43 * unit::deg;
}  // namespace space::consts

namespace space::consts {
  constexpr double G = 1;
  constexpr double C = 299792.458 * unit::kms;
}  // namespace space::consts

#endif
