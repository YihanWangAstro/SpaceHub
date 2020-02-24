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
 * @file orbits/orbits.hpp
 *
 * Header file.
 */
#ifndef SPACEHUB_ORBITS_HPP
#define SPACEHUB_ORBITS_HPP

#include <cmath>
#include <variant>
#include "../IO.hpp"
#include "../core-computation.hpp"
#include "../rand-generator.hpp"
#include "../vector/vector3.hpp"
/**
 * @namespace space::orbit
 *
 * Documentation
 */
namespace space::orbit {

/*---------------------------------------------------------------------------*\
     Class OrbitArgs Declaration
\*---------------------------------------------------------------------------*/
/**
 * @brief Enum of kepler orbit type. Possible value: Ellipse, Parabola, Hyperbola, None.
 */
enum class OrbitType { Ellipse, Parabola, Hyperbola, None };

/**
 * @brief A place holder that indicates one of the three angles in orbital parameters will be randomly generaged.
 */
struct RandomIndicator {
} isotherm;

/**
 * @brief Orbital parameters of the Kepler orbit.
 *@tparam Real Floating point like type.
 */
template <typename Real>
struct KeplerOrbit {
 public:
  /**
   * @brief Floating point like type.
   */
  using Scalar = Real;

  /**
   * @brief Variant that can hold both a Scalar OR a random indictor.
   */
  using Variant = std::variant<Real, RandomIndicator>;
  /**
   *  @brief Mass of the primary object.
   */
  Scalar m1;
  /**
   * @brief Mass of the secondary object.
   */
  Scalar m2;
  /**
   * @brief Semi-latus rectum of the orbit  a(1-e^2)  .
   *
   * We don't use semi-major axis for this general orbital type because the semi-major axis for parabolic orbit is
   * undefined.
   */
  Scalar p;
  /**
   *  @brief Eccentricity of the orbit.
   */
  Scalar e;
  /**
   *  @brief Orbit inclination.
   */
  Scalar i;
  /**
   *  @brief Longitude of the ascending node.
   */
  Scalar Omega;
  /**
   *  @brief Argument of periapsis.
   */
  Scalar omega;
  /**
   *  @brief True anomaly.
   */
  Scalar nu;
  /**
   *  @brief Orbit type.
   */
  OrbitType orbit_type;

  SPACEHUB_MAKE_CONSTRUCTORS(KeplerOrbit, default, default, default, default, default);

  /**
   * @brief Construct a new Orbit Args object from orbital parameters
   *
   * @param[in] m_1 Mass of the primary object.
   * @param[in] m_2 Mass of the secondary object.
   * @param[in] semi_latus_rectum Semi-latus rectum.
   * @param[in] eccentricity Eccentricity.
   * @param[in] inclination Inclination.
   * @param[in] longitude_of_ascending_node Longitude of the ascending node.
   * @param[in] argument_of_periapsis Argument of the periapsis.
   * @param[in] true_anomaly True anomaly.
   */
  KeplerOrbit(Scalar m_1, Scalar m_2, Scalar semi_latus_rectum, Scalar eccentricity, Variant inclination,
              Variant longitude_of_ascending_node, Variant argument_of_periapsis, Variant true_anomaly);

  /**
   * @brief Suffle the inclination.
   */
  inline void shuffle_i();

  /**
   * @brief Suffle the Longitude of the ascending node.
   */
  inline void shuffle_Omega();

  /**
   * @brief Suffle the Argument of periapsis.
   */
  inline void shuffle_omega();

  /**
   * @brief Suffle the true anomaly.
   */
  inline void shuffle_nu();

  /**
   * @brief Write the orbit to an ouput stream.
   *
   * @param[out] os Output stream
   * @param[in] obt Orbit parameters.
   * @return std::ostream& Output stream.
   */
  friend std::ostream &operator<<(std::ostream &os, KeplerOrbit const &obt) {
    space::display(os, obt.m1, obt.m2, obt.p, obt.e, obt.i, obt.Omega, obt.omega, obt.nu);
    return os;
  }
};

/**
 * @brief Alias of OrbitArgs<double>.
 */
using Kepler = KeplerOrbit<double>;

/*---------------------------------------------------------------------------*\
     Class HyperOrbit Declaration
\*---------------------------------------------------------------------------*/
/**
 * @brief Enum type that indicates the trajectory is hyperbolically incident in or hyperbolically eject out.
 */
enum class Hyper { in, out };

/**
 * @brief Derived class of Kepler orbit. Hyperbolic orbit.
 */
struct HyperOrbit : public KeplerOrbit<double> {
 public:
  /**
   * @brief Variant that can hold both a Scalar OR a random indictor.
   */
  using Variant = typename KeplerOrbit<double>::Variant;
  /**
   * @brief Floating point like type.
   */
  using Scalar = typename KeplerOrbit<double>::Scalar;

  SPACEHUB_MAKE_CONSTRUCTORS(HyperOrbit, delete, default, default, default, default);

  /**
   * @brief Construct a new Hyper Orbit object from scattering parameter b:impact parameter, v_inf:velocity at infinity.
   *
   * @param[in] m_1 Mass of the primary object: Stayed object.
   * @param[in] m_2 Mass of the secondary object: incident object.
   * @param[in]  v_inf Velocity at infinity.
   * @param[in]  b Impact parameter.
   * @param[in]  inclination Orbit inclination
   * @param[in]  longitude_of_ascending_node Longitude of the ascending node.
   * @param[in]  argument_of_periapsis Argument of the periapsis.
   * @param[in]  r Distance between the two objects.
   * @param[in]  in_out  Indicator of incident in or ejected out
   */
  HyperOrbit(Scalar m_1, Scalar m_2, Scalar v_inf, Scalar b, Variant inclination, Variant longitude_of_ascending_node,
             Variant argument_of_periapsis, Scalar r, Hyper in_out = Hyper::in);

  /**
   * @brief impact parameter.
   */
  Scalar b{0};
};

/*---------------------------------------------------------------------------*\
     Class EllipOrbit Declaration
\*---------------------------------------------------------------------------*/
/**
 * @brief Derived class of Kepler orbit. Elliptical orbit.
 */
struct EllipOrbit : public KeplerOrbit<double> {
 public:
  /**
   * @brief Variant that can hold both a Scalar OR a random indictor.
   */
  using Variant = typename KeplerOrbit<double>::Variant;
  /**
   * @brief Floating point like type.
   */
  using Scalar = typename KeplerOrbit<double>::Scalar;

  SPACEHUB_MAKE_CONSTRUCTORS(EllipOrbit, delete, default, default, default, default);

  /**
   * @brief Construct a new Elliptical Orbit object from orbital parameters.
   *
   * @param[in] m_1 Mass of the primary object.
   * @param[in] m_2 Mass of the secondary object.
   * @param[in] semi_major_axis Semi-major axis.
   * @param[in] eccentricity Eccentricity.
   * @param[in] inclination Inclination.
   * @param[in] longitude_of_ascending_node Longitude of the ascending node.
   * @param[in] argument_of_periapsis Argument of periapsis.
   * @param[in] true_anomaly True anomaly.
   */
  EllipOrbit(Scalar m_1, Scalar m_2, Scalar semi_major_axis, Scalar eccentricity, Variant inclination,
             Variant longitude_of_ascending_node, Variant argument_of_periapsis, Variant true_anomaly);

  /**
   * @brief Semi-major axis.
   */
  Scalar a{0};
};

/*---------------------------------------------------------------------------*\
     Help functions
\*---------------------------------------------------------------------------*/
template <typename Scalar>
Scalar myacos(Scalar x) {
  return acos(math::in_range(-1.0, x, 1.0));
}

template <typename Scalar>
constexpr auto semi_latus_rectum(Scalar a, Scalar e) {
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

/**
 * @brief Calculate the corresponding true anomaly of the eccentric anomaly.
 *
 * @tparam Scalar Floating point like type.
 * @param E_anomaly Eccentric anomaly.
 * @param e Eccentricity.
 * @return Scalar True anomaly.
 */
template <typename Scalar>
Scalar E_anomaly_to_T_anomaly(Scalar E_anomaly, Scalar e) {
  if (0 <= e && e < 1)
    return 2 * atan2(sqrt(1 + e) * sin(E_anomaly * 0.5), sqrt(1 - e) * cos(0.5 * E_anomaly));
  else if (e > 1)
    return 2 * atan2(sqrt(1 + e) * sinh(E_anomaly * 0.5), sqrt(e - 1) * cosh(0.5 * E_anomaly));
  else if (math::iseq(e, 1.0))
    return 2 * atan(0.5 * E_anomaly);
  else {
    spacehub_abort("Eccentricity cannot be negative, Nan or inf!");
  }
}

/**
 * @brief Calculate the corresponding eccentric anomaly of the mean anomaly.
 *
 * @tparam Scalar Floating point like type.
 * @param M_anomaly Mean anomaly.
 * @param e Eccentricity.
 * @return Scalar Eccentric anomaly.
 */
template <typename Scalar>
Scalar M_anomaly_to_E_anomaly(Scalar M_anomaly, Scalar e) {
  if (fabs(M_anomaly) <= math::epsilon<Scalar>::value) {
    return 0;
  }

  if (0 <= e && e < 1)
    return math::root_bisection([=](Scalar x) -> Scalar { return (x - e * sin(x) - M_anomaly) / (1 - e * cos(x)); },
                                -space::consts::pi, space::consts::pi);
  else if (e > 1)
    return math::root_bisection([=](Scalar x) -> Scalar { return (e * sinh(x) - x - M_anomaly) / (e * cosh(x) - 1); },
                                -space::consts::pi, space::consts::pi);
  else if (fabs(e - 1) < math::epsilon<Scalar>::value)
    return math::root_bisection([=](Scalar x) -> Scalar { return (x + x * x * x / 3 - M_anomaly) / (1 + x * x); },
                                -space::consts::pi, space::consts::pi);
  else {
    spacehub_abort("Eccentricity cannot be negative, Nan or inf!");
  }
}

/**
 * @brief Calculate the corresponding true anomaly of the mean anomaly.
 *
 * @tparam Scalar Floating point like type.
 * @param M_anomaly Mean anomaly.
 * @param e Eccentricity.
 * @return Scalar True anomaly.
 */
template <typename Scalar>
Scalar M_anomaly_to_T_anomaly(Scalar M_anomaly, Scalar e) {
  return E_anomaly_to_T_anomaly(M_anomaly_to_E_anomaly(M_anomaly, e), e);
}

/**
 * @brief Calculate the corresponding eccentric anomaly of the true anomaly.
 *
 * @tparam Scalar Floating point like type.
 * @param T_anomaly True anomaly.
 * @param e Eccentricity.
 * @return Scalar Eccentric anomaly.
 */
template <typename Scalar>
Scalar T_anomaly_to_E_anomaly(Scalar T_anomaly, Scalar e) {
  if (math::iseq(e, 1.0)) {
    return tan(0.5 * T_anomaly);
  } else if (0 <= e && e < 1) {
    auto cos_T = cos(T_anomaly);
    return acos((e + cos_T) / (1 + e * cos_T));
  } else if (e > 1) {
    auto cos_T = cos(T_anomaly);
    return acosh((e + cos_T) / (1 + e * cos_T));
  } else {
    spacehub_abort("Eccentricity cannot be negative, Nan or inf!");
  }
}

/**
 * @brief Calculate the corresponding mean anomaly of the eccentric anomaly.
 *
 * @tparam Scalar Floating point like type.
 * @param E_anomaly Eccentric anomaly.
 * @param e Eccentricity.
 * @return Scalar Mean anomaly.
 */
template <typename Scalar>
Scalar E_anomaly_to_M_anomaly(Scalar E_anomaly, Scalar e) {
  if (0 <= e && e < 1) {
    return E_anomaly - e * sin(E_anomaly);
  } else if (e > 1) {
    return e * sinh(E_anomaly) - E_anomaly;
  } else if (math::iseq(e, 1.0)) {
    return E_anomaly + E_anomaly * E_anomaly * E_anomaly / 3.0;
  } else {
    spacehub_abort("Eccentricity cannot be negative, Nan or inf!");
  }
}

/**
 * @brief Calculate the corresponding mean anomaly of the true anomaly.
 *
 * @tparam Scalar Floating point like type.
 * @param T_anomaly True anomaly.
 * @param e Eccentricity.
 * @return Scalar Mean anomaly.
 */
template <typename Scalar>
Scalar T_anomaly_to_M_anomaly(Scalar T_anomaly, Scalar e) {
  E_anomaly_to_M_anomaly(T_anomaly_to_E_anomaly(T_anomaly, e), e);
}
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

/*---------------------------------------------------------------------------*\
     Class OrbitArgs Implementation
\*---------------------------------------------------------------------------*/
template <typename Real>
KeplerOrbit<Real>::KeplerOrbit(Scalar m_1, Scalar m_2, Scalar semi_latus_rectum, Scalar eccentricity,
                               Variant inclination, Variant longitude_of_ascending_node, Variant argument_of_periapsis,
                               Variant true_anomaly) {
  if (semi_latus_rectum < 0) spacehub_abort("Semi-latus rectum cannot be negative");

  orbit_type = classify_orbit(eccentricity);

  if (orbit_type == OrbitType::None) {
    spacehub_abort("Eccentricity cannot be negative or NaN!");
  }

  m1 = m_1;
  m2 = m_2;
  p = semi_latus_rectum;
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
void KeplerOrbit<Real>::shuffle_i() {
  i = acos(random::Uniform(-1, 1));
}

template <typename Real>
void KeplerOrbit<Real>::shuffle_Omega() {
  Omega = random::Uniform(-consts::pi, consts::pi);
}

template <typename Real>
void KeplerOrbit<Real>::shuffle_omega() {
  omega = random::Uniform(-consts::pi, consts::pi);
}

template <typename Real>
void KeplerOrbit<Real>::shuffle_nu() {
  if (orbit_type == OrbitType::Ellipse) {
    Scalar M = random::Uniform(-consts::pi, consts::pi);
    Scalar E = orbit::M_anomaly_to_E_anomaly(M, e);
    nu = orbit::E_anomaly_to_T_anomaly(E, e);
  } else {
    spacehub_abort("Only elliptical orbit provides random anomaly method at this moment!");
  }
}

/*---------------------------------------------------------------------------*\
     Class HyperOrbit Implementation
\*---------------------------------------------------------------------------*/
HyperOrbit::HyperOrbit(Scalar m_1, Scalar m_2, Scalar v_inf, Scalar b, Variant inclination,
                       Variant longitude_of_ascending_node, Variant argument_of_periapsis, Scalar r, Hyper in_out)
    : KeplerOrbit<double>(m_1, m_2, 0, 0, inclination, longitude_of_ascending_node, argument_of_periapsis, 0) {
  this->orbit_type = OrbitType::Hyperbola;
  Scalar u = space::consts::G * (m_1 + m_2);
  Scalar a = -u / (v_inf * v_inf);
  this->e = sqrt(1 + b * b / (a * a));
  this->p = a * (1 - e * e);
  this->nu = -acos((p - r) / (e * r));
  this->b = b;
  if (in_out == Hyper::out) {
    this->nu *= -1;
  }
}

/*---------------------------------------------------------------------------*\
     Class EllipOrbit Implementation
\*---------------------------------------------------------------------------*/
EllipOrbit::EllipOrbit(Scalar m_1, Scalar m_2, Scalar semi_major_axis, Scalar eccentricity, Variant inclination,
                       Variant longitude_of_ascending_node, Variant argument_of_periapsis, Variant true_anomaly)
    : KeplerOrbit<double>(m_1, m_2, semi_major_axis * (1 - eccentricity * eccentricity), eccentricity, inclination,
                          longitude_of_ascending_node, argument_of_periapsis, true_anomaly) {
  if (this->orbit_type != OrbitType::Ellipse) {
    spacehub_abort("The given parameters don't give an elliptic orbit.");
  }
  a = semi_major_axis;
}

/*---------------------------------------------------------------------------*\
    Help functions
\*---------------------------------------------------------------------------*/
template <typename Vector, typename Scalar>
void orbit_to_coord(const KeplerOrbit<Scalar> &args, Vector &pos, Vector &vel) {
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

template <typename Vector, typename Scalar>
void coord_to_orbit(Scalar m1, Scalar m2, const Vector &dr, const Vector &dv, KeplerOrbit<Scalar> &args) {
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

/**
 * @brief Tranfer relative position and velocity between two particles to Kepler orbit parameters.
 *
 * @tparam Vector
 * @tparam Scalar
 * @param m1
 * @param m2
 * @param dr
 * @param dv
 * @return auto
 */
template <typename Vector, typename Scalar>
auto coord_to_orbit(Scalar m1, Scalar m2, const Vector &dr, const Vector &dv) {
  Vector L = cross(dr, dv);
  Vector N = cross(Vector(0, 0, 1.0), L);
  Scalar r = norm(dr);
  Scalar n = norm(N);
  Scalar l = norm(L);
  Scalar rv = dot(dr, dv);
  Scalar u = (m1 + m2) * consts::G;
  Vector E = (dr * (norm2(dv) - u * re_norm(dr)) - dv * rv) / u;

  KeplerOrbit<Scalar> args;

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
  return args;
}

/**
 * @brief Transfer Kepler orbit parameters to relative position and velocity between two component in orbit.
 *
 * @tparam Scalar
 * @param args
 * @return auto
 */
template <typename Scalar>
auto orbit_to_coord(KeplerOrbit<Scalar> const &args) {
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
inline auto calc_runge_lenz_vector(Scalar u, Vector const &dr, Vector const &dv) {
  return (dr * (norm2(dv) - u * re_norm(dr)) - dv * dot(dr, dv)) / u;
}

template <typename Scalar>
inline auto calc_runge_lenz_vector(Scalar u, Scalar dx, Scalar dy, Scalar dz, Scalar dvx, Scalar dvy, Scalar dvz) {
  using Vector = Vec3<Scalar>;
  return calc_runge_lenz_vector(u, Vector(dx, dy, dz), Vector(dvx, dvy, dvz));
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

template <typename Vector, typename Scalar>
auto calc_a_RL_vector(Scalar u, Vector const &dr, Vector const &dv) {
  Scalar r = norm(dr);
  Scalar v = norm(dv);
  Scalar v2 = v * v;
  Scalar vr = dot(dr, dv);
  Scalar vdfs = v2 - u / r;
  Scalar a = -u / (v2 - 2 * u / r);
  Vector e =((dr * vdfs - dv * vr) / u);
  return std::make_tuple(a, e);
}

template <typename Scalar>
inline auto calc_a_RL_vector(Scalar u, Scalar dx, Scalar dy, Scalar dz, Scalar dvx, Scalar dvy, Scalar dvz) {
  using Vector = Vec3<Scalar>;
  return calc_a_RL_vector(u, Vector(dx, dy, dz), Vector(dvx, dvy, dvz));
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
inline auto period(KeplerOrbit<Scalar> const &args) {
  return period(args.m1, args.m2, args.p / (1 - args.e * args.e));
}

template <typename Scalar>
auto time_to_periapsis(OrbitType obt_type, Scalar u, Scalar a, Scalar M_anomaly) {
  if (obt_type == OrbitType::Ellipse) {
    return sqrt(a * a * a / u) * M_anomaly;
  } else if (obt_type == OrbitType::Parabola) {
    return 0.5 * sqrt(a * a * a / u) * M_anomaly;
  } else if (obt_type == OrbitType::Hyperbola) {
    return sqrt(-a * a * a / u) * M_anomaly;
  }
}

template <typename Scalar>
auto time_to_periapsis(KeplerOrbit<Scalar> const &args) {
  auto M = E_anomaly_to_M_anomaly(T_anomaly_to_E_anomaly(args.nu));
  auto u = consts::G * (args.m1 + args.m2);
  if (args.orbit_type == OrbitType::Ellipse) {
    auto a = args.p / (1 - args.e * args.e);
    return sqrt(a * a * a / u) * M;
  } else if (args.orbit_type == OrbitType::Parabola) {
    return 0.5 * sqrt(args.p * args.p * args.p / u) * M;
  } else if (args.orbit_type == OrbitType::Hyperbola) {
    auto a = args.p / (1 - args.e * args.e);
    return sqrt(-a * a * a / u) * M;
  }
}

template <typename Scalar>
auto tidal_factor(Scalar r, Scalar m_tot1, Scalar m_tot2, Scalar R1, Scalar R2) {
  auto ratio1 = R1 / r;
  auto ratio2 = R2 / r;

  return std::make_tuple(m_tot2 / m_tot1 * ratio1 * ratio1 * ratio1, m_tot1 / m_tot2 * ratio2 * ratio2 * ratio2);
}

template <typename Scalar>
auto tidal_radius(Scalar tidal_factor, Scalar m_tot1, Scalar m_tot2, Scalar R2) {
  return pow(m_tot1 / (tidal_factor * m_tot2), 1.0 / 3) * R2;
}

template <typename Scalar>
auto tidal_radius(Scalar tidal_factor, Scalar m1, Scalar m2, Scalar R1, Scalar R2) {
  auto r1 = tidal_radius(tidal_factor, m1, m2, R2);
  auto r2 = tidal_radius(tidal_factor, m2, m1, R1);
  return std::max(r1, r2);
}

}  // namespace space::orbit
#endif
