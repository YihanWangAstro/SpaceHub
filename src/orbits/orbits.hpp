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

/**
 *
 * @tparam Real
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

struct HyperOrbit : public OrbitArgs<double> {
 private:
  using Variant = std::variant<double, RandomIndicator>;

 public:
  using Scalar = double;

  HyperOrbit() = default;

  HyperOrbit(Scalar m_1, Scalar m_2, Scalar v_inf, Scalar b, Scalar r, Variant inclination,
             Variant longitude_of_ascending_node, Variant argument_of_periapsis, Hyper in_out = Hyper::in);
};

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

EllipOrbit::EllipOrbit(Scalar m_1, Scalar m_2, Scalar semi_major_axis, Scalar eccentricity, Variant inclination,
                       Variant longitude_of_ascending_node, Variant argument_of_periapsis, Variant true_anomaly)
    : OrbitArgs<double>(m_1, m_2, semi_major_axis * (1 - eccentricity * eccentricity), eccentricity, inclination,
                        longitude_of_ascending_node, argument_of_periapsis, true_anomaly) {
  if (this->orbit_type != OrbitType::Ellipse) {
    spacehub_abort("The given parameters don't give an elliptic orbit.");
  }
  a = semi_major_axis;
}

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

template <typename Particle, typename... Args>
auto cluster(Particle const &ptc1, Particle const &ptc2, Args const &... ptcs) {
  static_assert(calc::all(std::is_same_v<Args, Particle>...), "Type of the arguments must be same!");
  return std::initializer_list<Particle>{ptc1, ptc2, ptcs...};
}

template <typename Particle, typename... Args>
inline auto M_tot(Particle &&ptc, Args &&... args) {
  if constexpr (sizeof...(Args) != 0) {
    static_assert(calc::all(std::is_same_v<Args, Particle>...),
                  "Type of the 1st argument and the rest should be same!");
    return (args.mass + ... + ptc.mass);
  } else if constexpr (is_container_v<Particle>) {
    typename Particle::Scaler tot_m = 0;
    for (auto &p : ptc) {
      tot_m += p.mass;
    }
    return tot_m;
  } else {
    return ptc.mass;
  }
}

template <typename Particle, typename... Args>
inline auto COM_p(Particle const &ptc, Args const &... ptcs) {
  if constexpr (sizeof...(Args) != 0) {
    static_assert(calc::all(std::is_same_v<Args, Particle>...),
                  "Type of the 1st argument and the rest should be same!");
    auto tot_mass = (ptcs.mass + ... + ptc.mass);
    auto cm_pos = ((ptcs.mass * ptcs.pos) + ... + (ptc.mass * ptc.pos)) / tot_mass;

    return cm_pos;
  } else if constexpr (is_container_v<Particle>) {
    using sParticle = typename Particle::value_type;
    using Scalar = typename sParticle::Scalar;
    using Vector = typename sParticle::Vector;

    Scalar tot_mass{0};
    Vector cm_pos{0};

    for (auto const &p : ptc) {
      tot_mass += p.mass;
      cm_pos += p.mass * p.pos;
    }

    cm_pos /= tot_mass;

    return cm_pos;
  } else {
    return ptc.pos;
  }
}

template <typename Particle, typename... Args>
inline auto COM_v(Particle const &ptc, Args const &... ptcs) {
  if constexpr (sizeof...(Args) != 0) {
    static_assert(calc::all(std::is_same_v<Args, Particle>...),
                  "Type of the 1st argument and the rest should be same!");
    auto tot_mass = (ptcs.mass + ... + ptc.mass);
    auto cm_vel = ((ptcs.mass * ptcs.vel) + ... + (ptc.mass * ptc.vel)) / tot_mass;

    return cm_vel;
  } else if constexpr (is_container_v<Particle>) {
    using sParticle = typename Particle::value_type;
    using Scalar = typename sParticle::Scalar;
    using Vector = typename sParticle::Vector;

    Scalar tot_mass{0};
    Vector cm_vel{0};

    for (auto const &p : ptc) {
      tot_mass += p.mass;
      cm_vel += p.mass * p.vel;
    }

    cm_vel /= tot_mass;

    return cm_vel;
  } else {
    return ptc.vel;
  }
}

template <typename Vector, typename Particle, typename... Args>
void move_particles_pos(Vector const &centre_mass_pos, Particle &ptc, Args &... ptcs) {
  auto dp = centre_mass_pos - COM_p(ptc, ptcs...);
  if constexpr (sizeof...(Args) != 0) {
    static_assert(calc::all(std::is_same_v<Args, Particle>...),
                  "Type of the 1st argument and the rest should be same!");
    ((ptc.pos += dp), ..., (ptcs.pos += dp));
  } else if constexpr (is_container_v<Particle>) {
    for (auto &p : ptc) {
      p.pos += dp;
    }
  } else {
    ptc.pos += dp;
  }
}

template <typename Vector, typename Particle, typename... Args>
void move_particles_vel(Vector const &centre_mass_vel, Particle &ptc, Args &... ptcs) {
  auto dv = centre_mass_vel - COM_v(ptc, ptcs...);
  if constexpr (sizeof...(Args) != 0) {
    static_assert(calc::all(std::is_same_v<Args, Particle>...),
                  "Type of the 1st argument and the rest should be same!");
    ((ptc.vel += dv), ..., (ptcs.vel += dv));
  } else if constexpr (is_container_v<Particle>) {
    for (auto &p : ptc) {
      p.vel += dv;
    }
  } else {
    ptc.vel += dv;
  }
}

template <typename Vector, typename Particle, typename... Args>
void move_particles(Vector const &centre_mass_pos, Vector const &centre_mass_vel, Particle &ptc, Args &... ptcs) {
  auto dp = centre_mass_pos - COM_p(ptc, ptcs...);
  auto dv = centre_mass_vel - COM_v(ptc, ptcs...);
  if constexpr (sizeof...(Args) != 0) {
    static_assert(calc::all(std::is_same_v<Args, Particle>...),
                  "Type of the 1st argument and the rest should be same!");
    ((ptc.pos += dp, ptc.vel += dv), ..., (ptcs.pos += dp, ptcs.vel += dv));
  } else if constexpr (is_container_v<Particle>) {
    for (auto &p : ptc) {
      p.pos += dp;
      p.vel += dv;
    }
  } else {
    ptc.pos += dp;
    ptc.vel += dv;
  }
}

template <typename Scalar, typename Particle, typename... Args>
void move_particles(OrbitArgs<Scalar> const &args, Particle &ptc, Args &... ptcs) {
  static_assert(calc::all(std::is_same_v<Args, Particle>...), "Type of the 2nd argument and the rest should be same!");
  auto [cm_pos, cm_vel] = orbit_to_coord(args);
  move_particles(cm_pos, cm_vel, ptc, ptcs...);
}

template <typename Particle, typename... Args>
void move_to_COM_frame(Particle &ptc, Args &... ptcs) {
  using Vector = decltype(COM_v(ptc, ptcs));
  move_particles(Vector{0, 0, 0}, Vector{0, 0, 0}, ptc, ptcs...);
}

template <typename Scalar>
inline constexpr auto semi_latus_rectum(Scalar a, Scalar e) {
  return a * (1 - e * e);
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

template <typename Particle1, typename Particle2>
inline auto calc_eccentricity(Particle1 const &p1, Particle2 const &p2) {
  auto m1 = M_tot(p1);
  auto m2 = M_tot(p2);
  auto dp = COM_p(p1) - COM_p(p2);
  auto dv = COM_v(p1) - COM_v(p2);
  return calc_eccentricity(consts::G * (m1 + m2), dp, dv);
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

template <typename Particle1, typename Particle2>
inline auto calc_semi_major_axis(Particle1 const &p1, Particle2 const &p2) {
  auto m1 = M_tot(p1);
  auto m2 = M_tot(p2);
  auto dp = COM_p(p1) - COM_p(p2);
  auto dv = COM_v(p1) - COM_v(p2);
  return calc_semi_major_axis(consts::G * (m1 + m2), dp, dv);
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

template <typename Particle>
inline auto calc_a_e(Particle const &p1, Particle const &p2) {
  auto m1 = M_tot(p1);
  auto m2 = M_tot(p2);
  auto dp = COM_p(p1) - COM_p(p2);
  auto dv = COM_v(p1) - COM_v(p2);
  return calc_a_e(consts::G * (m1 + m2), dp, dv);
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

template <typename Particle>
inline auto period(Particle const &p1, Particle const &p2) {
  auto m1 = M_tot(p1);
  auto m2 = M_tot(p2);
  return period(m1, m2, calc_semi_major_axis(p1, p2));
}

template <typename Particle, typename... Args>
auto E_k(Particle &&ptc, Args &&... args) {
  using Scalar = typename Particle::Scalar;
  using Vector = typename Particle::Vector;
  Scalar kinetic_energy = 0;

  auto com_vel = COM_v(std::forward<Particle>(ptc), std::forward<Args>(args)...);
  move_particles_vel(Vector{0, 0, 0}, std::forward<Particle>(ptc), std::forward<Args>(args)...);

  if constexpr (sizeof...(Args) != 0) {
    static_assert(calc::all(std::is_same_v<Args, Particle>...),
                  "Type of the 1st argument and the rest should be same!");

    kinetic_energy = ((args.mass * norm2(args.vel)) + ... + (ptc.mass * norm2(ptc.vel)));
  } else if constexpr (is_container_v<Particle>) {
    for (auto &p : ptc) {
      kinetic_energy += p.mass * norm2(p.vel);
    }
  }

  move_particles_vel(com_vel, std::forward<Particle>(ptc), std::forward<Args>(args)...);
  return kinetic_energy;
}

template <typename Particle, typename... Args>
auto E_p(Particle &&ptc, Args &&... args) {
  if constexpr (sizeof...(Args) != 0) {
    static_assert(calc::all(std::is_same_v<Args, Particle>...),
                  "Type of the 1st argument and the rest should be same!");
    std::initializer_list<Particle> list{ptc, args...};
    return E_p(list);
  } else if constexpr (is_container_v<Particle>) {
    typename Particle::Scaler potential_energy = 0;

    for (auto i = ptc.begin(); i < ptc.end(); ++i) {
      for (auto j = i + 1; j < ptc.end(); ++j) {
        potential_energy -= consts::G * i->mass * j->mass / distance(i->pos, j->pos);
      }
    }
    return potential_energy;
  } else {
    return static_cast<typename Particle::Scaler>(0);
  }
}

CREATE_MEMBER_CHECK(mass);
CREATE_MEMBER_CHECK(radius);

template <typename Particle, typename... Args>
auto E_tot(Particle &&ptc, Args &&... args) {
  return E_k(std::forward<Particle>(ptc), std::forward<Args>(args)...) +
         E_p(std::forward<Particle>(ptc), std::forward<Args>(args)...);
}

template <typename Particle>
inline auto M_rdc(Particle &&m1, Particle &&m2) {
  auto tot_mass1 = M_tot(m1);
  auto tot_mass2 = M_tot(m2);
  return tot_mass1 * tot_mass2 / (tot_mass1 + tot_mass2);
}

template <typename Particle, typename... Args>
auto M_cluster_rdc(Particle &&ptc, Args &&... args) {
  if constexpr (sizeof...(Args) != 0) {
    static_assert(calc::all(std::is_same_v<Args, Particle>...),
                  "Type of the 1st argument and the rest should be same!");
    return M_cluster_rdc(cluster(ptc, args...));
  } else if constexpr (is_container_v<Particle>) {
    typename Particle::Scaler mass_product = 0;

    size_t particle_num = ptc.size();

    if (particle_num > 1) {
      for (auto i = ptc.begin(); i < ptc.end(); ++i) {
        for (auto j = i + 1; j < ptc.end(); ++j) {
          mass_product += i->mass * j->mass;
        }
      }
      return mass_product / M_tot(ptc);
    } else if (particle_num == 1) {
      return ptc.begin()->mass;
    } else {
      return static_cast<typename Particle::Scaler>(0);
    }
  } else {
    return ptc.mass;
  }
}

template <typename Particle, typename... Args>
auto cluster_size(Particle &&ptc, Args &&... args) {
  if constexpr (sizeof...(Args) != 0) {
    static_assert(calc::all(std::is_same_v<Args, Particle>...),
                  "Type of the 1st argument and the rest should be same!");
    return cluster_size(cluster(ptc, args...));
  } else if constexpr (is_container_v<Particle>) {
    typename Particle::Scaler R_max = 0;
    size_t particle_num = ptc.size();

    if (particle_num > 1) {
      for (auto i = ptc.begin(); i < ptc.end(); ++i) {
        for (auto j = i + 1; j < ptc.end(); ++j) {
          auto r = distance(i->pos, j->pos);
          if (r > R_max) {
            R_max = r;
          }
        }
      }
      return R_max;
    } else if (particle_num == 1) {
      return cluster_size(*(ptc.begin()));
    } else {
      return static_cast<typename Particle::value_type::Scaler>(0);
    }
  } else {
    if constexpr (HAS_MEMBER(Particle, radius)) {
      return ptc.radius;
    } else {
      return static_cast<typename Particle::Scaler>(0);
    }
  }
}

template <typename T1, typename T2, typename Scalar>
auto E_tid(T1 &&m1, T2 &&m2, Scalar R2) {
  auto r = distance(orbit::COM_p(m1), orbit::COM_p(m2));
  auto m_tot1 = orbit::M_tot(m1);
  auto m_tot2 = orbit::M_tot(m2);

  return -consts::G * m_tot1 * m_tot2 * R2 / (r * r);
}

template <typename T1, typename T2>
auto E_tid(T1 &&m1, T2 &&m2) {
  auto R2 = orbit::cluster_size(m2);
  return E_tid(std::forward<T1>(m1), std::forward<T2>(m2), R2);
}

template <typename T1, typename T2, typename Scalar>
auto tidal_factor(T1 &&m1, T2 &&m2, Scalar R2) {
  auto r = distance(orbit::COM_p(m1), orbit::COM_p(m2));
  auto m_tot1 = orbit::M_tot(m1);
  auto m_tot2 = orbit::M_tot(m2);

  auto ratio = R2 / r;

  return m_tot1 / m_tot2 * ratio * ratio * ratio;
}

template <typename T1, typename T2>
auto tidal_factor(T1 &&m1, T2 &&m2) {
  auto R2 = orbit::cluster_size(m2);
  return tidal_factor(std::forward<T1>(m1), std::forward<T2>(m2), R2);
}

template <typename T1, typename T2, typename Scalar>
auto tidal_radius(Scalar tidal_factor, T1 &&m1, T2 &&m2, Scalar R2) {
  auto m_tot1 = orbit::M_tot(m1);
  auto m_tot2 = orbit::M_tot(m2);

  return pow(m_tot1 / (tidal_factor * m_tot2), 1.0 / 3) * R2;
}

template <typename T1, typename T2, typename Scalar>
auto tidal_radius(Scalar tidal_factor, T1 &&m1, T2 &&m2) {
  auto R2 = orbit::cluster_size(m2);
  return tidal_radius(tidal_factor, std::forward<T1>(m1), std::forward<T2>(m2), R2);
}

}  // namespace space::orbit
#endif
