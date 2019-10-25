/*---------------------------------------------------------------------------*\
        .-''''-.         |
       /        \        |
      /_        _\       |  SpaceHub: The Open Source N-body Toolkit
     // \  <>  / \\      |
     |\__\    /__/|      |  Website:  https://yihanwangastro.github.io/SpaceHub/
      \    ||    /       |
        \  __  /         |  Copyright (C) 2019 Yihan Wang
         '.__.'          |
---------------------------------------------------------------------
License
    This file is part of SpaceHub.
    SpaceHub is free software: you can redistribute it and/or modify it under
    the terms of the MIT License. SpaceHub is distributed in the hope that it
    will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the MIT License
    for more details. You should have received a copy of the MIT License along
    with SpaceHub.
\*---------------------------------------------------------------------------*/
/**
 * @file orbits.hpp
 *
 * Header file.
 */
#ifndef SPACEHUB_ORBITS_HPP
#define SPACEHUB_ORBITS_HPP

#include <cmath>
#include <variant>
#include "../core-computation.hpp"
#include "../rand-generator.hpp"
#include "../vector/vector3.hpp"
/**
 * @namespace space::orbit
 * Documentation for space
 */
namespace space::orbit {

template <typename Scalar>
inline Scalar myacos(Scalar x) {
  return acos(math::in_range(-1.0, x, 1.0));
}

template <typename Scalar>
inline constexpr auto semi_latus_rectum(Scalar a, Scalar e) {
  return a * (1 - e * e);
}

template <typename Vector, typename Scalar>
void euler_rotate(Vector &v, const Scalar phi, const Scalar theta, const Scalar psi) {
  Scalar sin_phi = sin(phi);
  Scalar cos_phi = cos(phi);
  Scalar sin_psi = sin(psi);
  Scalar cos_psi = cos(psi);
  Scalar sin_theta = sin(theta);
  Scalar cos_theta = cos(theta);

  Scalar x = v.x * (cos_phi * cos_psi - sin_phi * cos_theta * sin_psi) -
             v.y * (cos_phi * sin_psi + sin_phi * cos_theta * cos_psi) + v.z * (sin_phi * sin_theta);

  Scalar y = v.x * (sin_phi * cos_psi + cos_phi * cos_theta * sin_psi) -
             v.y * (sin_phi * sin_psi - cos_phi * cos_theta * cos_psi) - v.z * (cos_phi * sin_theta);

  Scalar z = v.x * sin_theta * sin_psi + v.y * sin_theta * cos_psi + v.z * cos_theta;

  v.x = x, v.y = y, v.z = z;
}

template <typename Scalar>
Scalar calc_true_anomaly(Scalar eccentric_anomaly, Scalar e) {
  if (0 <= e && e < 1)
    return 2 * atan2(sqrt(1 + e) * sin(eccentric_anomaly * 0.5), sqrt(1 - e) * cos(0.5 * eccentric_anomaly));
  else if (e > 1)
    return 2 * atan2(sqrt(1 + e) * sinh(eccentric_anomaly * 0.5), sqrt(e - 1) * cosh(0.5 * eccentric_anomaly));
  else if (math::iseq(e, 1.0))
    return 2 * atan(0.5 * eccentric_anomaly);
  else {
    spacehub_abort("Eccentricity cannot be negative, Nan or inf!");
  }
}

template <typename Scalar>
Scalar calc_eccentric_anomaly(Scalar mean_anomaly, Scalar e) {
  if (fabs(mean_anomaly) <= math::epsilon<Scalar>::value) {
    return 0;
  }

  if (0 <= e && e < 1)
    return math::root_bisection([=](Scalar x) -> Scalar { return (x - e * sin(x) - mean_anomaly) / (1 - e * cos(x)); },
                                -space::consts::pi, space::consts::pi);  // find this function in ownMath.h
  else if (e > 1)
    return math::root_bisection(
        [=](Scalar x) -> Scalar { return (e * sinh(x) - x - mean_anomaly) / (e * cosh(x) - 1); }, -space::consts::pi,
        space::consts::pi);
  else if (fabs(e - 1) < math::epsilon<Scalar>::value)
    return math::root_bisection([=](Scalar x) -> Scalar { return (x + x * x * x / 3 - mean_anomaly) / (1 + x * x); },
                                -space::consts::pi, space::consts::pi);
  else {
    spacehub_abort("Eccentricity cannot be negative, Nan or inf!");
  }
}

enum class OrbitType { Ellipse, Parabola, Hyperbola, None };

template <typename T>
constexpr OrbitType classify_orbit(T eccentricity) {
  if (0 <= eccentricity && eccentricity < 1) {
    return OrbitType::Ellipse;
  } else if (math::iseq(eccentricity, 1.0)) {
    return OrbitType::Parabola;
  } else if (eccentricity > 1) {
    return OrbitType::Hyperbola;
  } else {
    return OrbitType::None;
  }
}

struct RandomIndicator {
} isotherm;

/*---------------------------------------------------------------------------*\
     Class OrbitArgs Declaration
\*---------------------------------------------------------------------------*/
/**

 *@tparam Real
 */
template <typename Real>
struct OrbitArgs {
 private:
  using Variant = std::variant<Real, RandomIndicator>;

 public:
  using Scalar = Real;
  // Scalar u;//gravitational parameter
  Scalar m1;
  Scalar m2;
  Scalar e;      // eccentricity
  Scalar p;      // semi-latus rectum
  Scalar i;      // inclination
  Scalar Omega;  // longitude of the ascending node
  Scalar omega;  // argument of periapsis
  Scalar nu;     // true anomaly
  OrbitType orbit_type;

  OrbitArgs() = default;

  OrbitArgs(Scalar m_1, Scalar m_2, Scalar periastron, Scalar eccentricity, Variant inclination,
            Variant longitude_of_ascending_node, Variant argument_of_periapsis, Variant true_anomaly);

  inline void shuffle_i();

  inline void shuffle_Omega();

  inline void shuffle_omega();

  inline void shuffle_nu();

  friend std::ostream &operator<<(std::ostream &os, OrbitArgs const &obt) {
    space::display(os, obt.m1, obt.m2, obt.e, obt.p, obt.i, obt.Omega, obt.omega, obt.nu);
    return os;
  }
};

/*---------------------------------------------------------------------------*\
     Class OrbitArgs Implementation
\*---------------------------------------------------------------------------*/
template <typename Real>
OrbitArgs<Real>::OrbitArgs(Scalar m_1, Scalar m_2, Scalar periastron, Scalar eccentricity, Variant inclination,
                           Variant longitude_of_ascending_node, Variant argument_of_periapsis, Variant true_anomaly) {
  if (periastron < 0) spacehub_abort("Semi-latus rectum cannot be negative");

  orbit_type = classify_orbit(eccentricity);

  if (orbit_type == OrbitType::None) {
    spacehub_abort("Eccentricity cannot be negative or NaN!");
  }

  m1 = m_1;
  m2 = m_2;
  p = periastron;
  e = eccentricity;

  if (std::holds_alternative<Scalar>(inclination)) {
    i = std::get<Scalar>(inclination);
  } else {
    shuffle_i();
  }

  if (std::holds_alternative<Scalar>(longitude_of_ascending_node)) {
    Omega = std::get<Scalar>(longitude_of_ascending_node);
  } else {
    shuffle_Omega();
  }

  if (std::holds_alternative<Scalar>(argument_of_periapsis)) {
    omega = std::get<Scalar>(argument_of_periapsis);
  } else {
    shuffle_omega();
  }

  if (std::holds_alternative<Scalar>(true_anomaly)) {
    nu = std::get<Scalar>(true_anomaly);
  } else {
    shuffle_nu();
  }
}

template <typename Real>
void OrbitArgs<Real>::shuffle_i() {
  i = acos(random::Uniform(-1, 1));
}

template <typename Real>
void OrbitArgs<Real>::shuffle_Omega() {
  Omega = random::Uniform(-consts::pi, consts::pi);
}

template <typename Real>
void OrbitArgs<Real>::shuffle_omega() {
  omega = random::Uniform(-consts::pi, consts::pi);
}

template <typename Real>
void OrbitArgs<Real>::shuffle_nu() {
  if (orbit_type == OrbitType::Ellipse) {
    Scalar M = random::Uniform(-consts::pi, consts::pi);
    Scalar E = orbit::calc_eccentric_anomaly(M, e);
    nu = orbit::calc_true_anomaly(E, e);
  } else {
    spacehub_abort("Only elliptical orbit provides random anomaly method at this moment!");
  }
}

using Kepler = OrbitArgs<double>;

enum class Hyper { in, out };

/*---------------------------------------------------------------------------*\
     Class HyperOrbit Declaration
\*---------------------------------------------------------------------------*/
struct HyperOrbit : public OrbitArgs<double> {
 private:
  using Variant = std::variant<double, RandomIndicator>;

 public:
  using Scalar = double;

  HyperOrbit() = default;

  HyperOrbit(Scalar m_1, Scalar m_2, Scalar v_inf, Scalar b, Scalar r, Variant inclination,
             Variant longitude_of_ascending_node, Variant argument_of_periapsis, Hyper in_out = Hyper::in);
};

/*---------------------------------------------------------------------------*\
     Class HyperOrbit Implementaion
\*---------------------------------------------------------------------------*/
HyperOrbit::HyperOrbit(Scalar m_1, Scalar m_2, Scalar v_inf, Scalar b, Scalar r, Variant inclination,
                       Variant longitude_of_ascending_node, Variant argument_of_periapsis, Hyper in_out)
    : OrbitArgs<double>(m_1, m_2, 0, 0, inclination, longitude_of_ascending_node, argument_of_periapsis, 0) {
  this->orbit_type = OrbitType::Hyperbola;
  Scalar u = space::consts::G * (m_1 + m_2);
  Scalar a = -u / (v_inf * v_inf);
  this->e = sqrt(1 + b * b / (a * a));
  this->p = a * (1 - e * e);
  this->nu = -acos((p - r) / (e * r));
  if (in_out == Hyper::out) {
    this->nu *= -1;
  }
}

/*---------------------------------------------------------------------------*\
     Class EllipOrbit Declaration
\*---------------------------------------------------------------------------*/
struct EllipOrbit : public OrbitArgs<double> {
 private:
  using Variant = std::variant<double, RandomIndicator>;

 public:
  using Scalar = double;

  Scalar a{0};

  EllipOrbit() = default;

  EllipOrbit(Scalar m_1, Scalar m_2, Scalar semi_major_axis, Scalar eccentricity, Variant inclination,
             Variant longitude_of_ascending_node, Variant argument_of_periapsis, Variant true_anomaly);
};

/*---------------------------------------------------------------------------*\
     Class EllipOrbit Implementation
\*---------------------------------------------------------------------------*/
EllipOrbit::EllipOrbit(Scalar m_1, Scalar m_2, Scalar semi_major_axis, Scalar eccentricity, Variant inclination,
                       Variant longitude_of_ascending_node, Variant argument_of_periapsis, Variant true_anomaly)
    : OrbitArgs<double>(m_1, m_2, semi_major_axis * (1 - eccentricity * eccentricity), eccentricity, inclination,
                        longitude_of_ascending_node, argument_of_periapsis, true_anomaly) {
  if (this->orbit_type != OrbitType::Ellipse) {
    spacehub_abort("The given parameters don't give an elliptic orbit.");
  }
  a = semi_major_axis;
}

/*---------------------------------------------------------------------------*\
    Functions Implementation
\*---------------------------------------------------------------------------*/
template <typename Vector, typename Scalar>
void coord_to_orbit(Scalar m1, Scalar m2, const Vector &dr, const Vector &dv, OrbitArgs<Scalar> &args) {
  Vector L = cross(dr, dv);
  Vector N = cross(Vector(0, 0, 1.0), L);
  Scalar r = norm(dr);
  Scalar n = norm(N);
  Scalar l = norm(L);
  Scalar rv = dot(dr, dv);
  Scalar u = (m1 + m2) * consts::G;
  Vector E = (dr * (norm2(dv) - u * re_norm(dr)) - dv * rv) / u;

  args.m1 = m1;
  args.m2 = m2;
  args.e = norm(E);
  args.orbit_type = classify_orbit(args.e);

  if (args.orbit_type == OrbitType::Parabola) {
    Scalar a = -u * r / (r * norm2(dv) - 2.0 * u);
    args.p = semi_latus_rectum(a, args.e);
  } else {
    args.p = l * l / u;
  }

  args.i = acos(L.z / l);

  if (args.e != 0) {
    args.nu = math::sign(rv) * myacos(dot(E / args.e, dr / r));

    if (n != 0) {
      args.Omega = math::sign(N.y) * myacos(N.x / n);
      args.omega = math::sign(E.z) * myacos(dot(E / args.e, N / n));
    } else {
      args.omega = -math::sign(E.y) * myacos(-E.x / args.e);
      args.Omega = args.omega;
    }
  } else {
    if (n != 0) {
      args.Omega = math::sign(N.y) * myacos(N.x / n);
      args.omega = 0;
      Vector peri = cross(L, N);
      args.nu = -math::sign(dot(N, dr)) * myacos(dot(peri / norm(peri), dr / r));
    } else {
      args.Omega = args.omega = 0;
      args.nu = math::sign(dr.y) * acos(dot(Vector(1.0, 0, 0), dr / r));
    }
  }
}

template <typename Vector, typename Scalar>
void orbit_to_coord(OrbitArgs<Scalar> const &args, Vector &pos, Vector &vel) {
  Scalar u = (args.m1 + args.m2) * consts::G;

  Scalar sin_nu = sin(args.nu);
  Scalar cos_nu = cos(args.nu);

  Scalar r = args.p / (1 + args.e * cos_nu);
  Scalar v = sqrt(u / args.p);

  pos = r * Vector(cos_nu, sin_nu, 0);
  vel = v * Vector(-sin_nu, args.e + cos_nu, 0);

  orbit::euler_rotate(pos, args.Omega, args.i, args.omega + consts::pi);
  orbit::euler_rotate(vel, args.Omega, args.i, args.omega + consts::pi);
}

template <typename Scalar>
auto orbit_to_coord(OrbitArgs<Scalar> const &args) {
  using Vector = Vec3<Scalar>;

  Scalar u = (args.m1 + args.m2) * consts::G;

  Scalar sin_nu = sin(args.nu);
  Scalar cos_nu = cos(args.nu);

  Scalar r = args.p / (1 + args.e * cos_nu);
  Scalar v = sqrt(u / args.p);

  Vector pos = r * Vector(cos_nu, sin_nu, 0);
  Vector vel = v * Vector(-sin_nu, args.e + cos_nu, 0);

  orbit::euler_rotate(pos, args.Omega, args.i, args.omega + consts::pi);
  orbit::euler_rotate(vel, args.Omega, args.i, args.omega + consts::pi);

  return std::make_tuple(pos, vel);
}

template <typename Vector, typename Scalar>
inline Scalar calc_eccentricity(Scalar u, Vector const &dr, Vector const &dv) {
  return norm(dr * (norm2(dv) - u * re_norm(dr)) - dv * dot(dr, dv)) / u;
}

template <typename Scalar>
inline auto calc_eccentricity(Scalar u, Scalar dx, Scalar dy, Scalar dz, Scalar dvx, Scalar dvy, Scalar dvz) {
  using Vector = Vec3<Scalar>;
  return calc_eccentricity(u, Vector(dx, dy, dz), Vector(dvx, dvy, dvz));
}

template <typename Vector, typename Scalar>
inline Scalar calc_semi_major_axis(Scalar u, Vector const &dr, Vector const &dv) {
  Scalar r = norm(dr);
  return -u * r / (r * norm2(dv) - 2 * u);
}

template <typename Scalar>
inline Scalar calc_semi_major_axis(Scalar u, Scalar dx, Scalar dy, Scalar dz, Scalar dvx, Scalar dvy, Scalar dvz) {
  using Vector = Vec3<Scalar>;
  return calc_semi_major_axis(u, Vector(dx, dy, dz), Vector(dvx, dvy, dvz));
}

template <typename Vector, typename Scalar>
auto calc_a_e(Scalar u, Vector const &dr, Vector const &dv) {
  Scalar r = norm(dr);
  Scalar v = norm(dv);
  Scalar v2 = v * v;
  Scalar vr = dot(dr, dv);
  Scalar vdfs = v2 - u / r;
  Scalar a = -u / (v2 - 2 * u / r);
  Scalar e = norm((dr * vdfs - dv * vr) / u);
  return std::make_tuple(a, e);
}

template <typename Scalar>
inline auto calc_a_e(Scalar u, Scalar dx, Scalar dy, Scalar dz, Scalar dvx, Scalar dvy, Scalar dvz) {
  using Vector = Vec3<Scalar>;
  return calc_a_e(u, Vector(dx, dy, dz), Vector(dvx, dvy, dvz));
}

template <typename Scalar>
inline auto period(Scalar m1, Scalar m2, Scalar a) {
  if (a > 0) {
    return 2 * consts::pi * sqrt(a * a * a / ((m1 + m2) * consts::G));
  } else {
    spacehub_abort("Only elliptical orbit is periodic!");
  }
}

template <typename Scalar>
inline auto period(OrbitArgs<Scalar> const &args) {
  return period(args.m1, args.m2, args.p / (1 - args.e * args.e));
}

}  // namespace space::orbit
#endif
