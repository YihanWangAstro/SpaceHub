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
 * @file operations.hpp
 *
 * Header file.
 */
#ifndef SPACEHUB_OPERATIONS_HPP
#define SPACEHUB_OPERATIONS_HPP

#include "../dev-tools.hpp"
#include "orbits.hpp"

namespace space::orbit {

template <typename Particle, typename... Args>
auto cluster(Particle const &ptc1, Particle const &ptc2, Args const &... ptcs) {
  static_assert(calc::all(std::is_same_v<Args, Particle>...), "Type of the arguments must be same!");
  return std::initializer_list<Particle>{ptc1, ptc2, ptcs...};
}

template <typename Cluster>
inline auto M_tot(Cluster &&ptc) {
  if constexpr (is_container_v<Cluster>) {
    typename Cluster::value_type::Scalar tot_m = 0;
    for (auto &p : ptc) {
      tot_m += p.mass;
    }
    return tot_m;
  } else {
    return ptc.mass;
  }
}

template <typename Particle, typename... Args>
inline auto M_tot(Particle &&ptc, Args &&... args) {
  static_assert(calc::all(std::is_same_v<Args, Particle>...), "Type of the 1st argument and the rest should be same!");
  return (args.mass + ... + ptc.mass);
}

template <typename Cluster>
inline auto COM_p(Cluster const &ptc) {
  if constexpr (is_container_v<Cluster>) {
    using Particle = typename Cluster::value_type;
    using Scalar = typename Particle::Scalar;
    using Vector = typename Particle::Vector;
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
inline auto COM_p(Particle const &ptc, Args const &... ptcs) {
  static_assert(calc::all(std::is_same_v<Args, Particle>...), "Type of the 1st argument and the rest should be same!");
  decltype(ptc.mass) tot_mass = (ptcs.mass + ... + ptc.mass);
  decltype(ptc.pos) cm_pos = ((ptcs.mass * ptcs.pos) + ... + (ptc.mass * ptc.pos)) / tot_mass;
  return cm_pos;
}

template <typename Cluster>
inline auto COM_v(Cluster const &ptc) {
  if constexpr (is_container_v<Cluster>) {
    using Particle = typename Cluster::value_type;
    using Scalar = typename Particle::Scalar;
    using Vector = typename Particle::Vector;

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

template <typename Particle, typename... Args>
inline auto COM_v(Particle const &ptc, Args const &... ptcs) {
  static_assert(calc::all(std::is_same_v<Args, Particle>...), "Type of the 1st argument and the rest should be same!");
  decltype(ptc.mass) tot_mass = (ptcs.mass + ... + ptc.mass);
  decltype(ptc.vel) cm_vel = ((ptcs.mass * ptcs.vel) + ... + (ptc.mass * ptc.vel)) / tot_mass;

  return cm_vel;
}

template <typename Cluster1, typename Cluster2>
inline auto M_rdc(Cluster1 &&m1, Cluster2 &&m2) {
  auto tot_mass1 = M_tot(m1);
  auto tot_mass2 = M_tot(m2);
  return tot_mass1 * tot_mass2 / (tot_mass1 + tot_mass2);
}

template <typename Vector, typename Cluster>
void move_particles_pos(Vector const &centre_mass_pos, Cluster &ptc) {
  auto dp = centre_mass_pos - COM_p(ptc);
  if constexpr (is_container_v<Cluster>) {
    for (auto &p : ptc) {
      p.pos += dp;
    }
  } else {
    ptc.pos += dp;
  }
}

template <typename Vector, typename Particle, typename... Args>
void move_particles_pos(Vector const &centre_mass_pos, Particle &ptc, Args &... ptcs) {
  static_assert(calc::all(std::is_same_v<Args, Particle>...), "Type of the 1st argument and the rest should be same!");
  auto dp = centre_mass_pos - COM_p(ptc, ptcs...);
  ((ptc.pos += dp), ..., (ptcs.pos += dp));
}

template <typename Vector, typename Cluster>
void move_particles_vel(Vector const &centre_mass_vel, Cluster &ptc) {
  auto dv = centre_mass_vel - COM_v(ptc);
  if constexpr (is_container_v<Cluster>) {
    for (auto &p : ptc) {
      p.vel += dv;
    }
  } else {
    ptc.vel += dv;
  }
}

template <typename Vector, typename Particle, typename... Args>
void move_particles_vel(Vector const &centre_mass_vel, Particle &ptc, Args &... ptcs) {
  static_assert(calc::all(std::is_same_v<Args, Particle>...), "Type of the 1st argument and the rest should be same!");
  auto dv = centre_mass_vel - COM_v(ptc, ptcs...);
  ((ptc.vel += dv), ..., (ptcs.vel += dv));
}

template <typename Vector, typename Cluster>
void move_particles(Vector const &centre_mass_pos, Vector const &centre_mass_vel, Cluster &ptc) {
  auto dp = centre_mass_pos - COM_p(ptc);
  auto dv = centre_mass_vel - COM_v(ptc);
  if constexpr (is_container_v<Cluster>) {
    for (auto &p : ptc) {
      p.pos += dp;
      p.vel += dv;
    }
  } else {
    ptc.pos += dp;
    ptc.vel += dv;
  }
}

template <typename Vector, typename Particle, typename... Args>
void move_particles(Vector const &centre_mass_pos, Vector const &centre_mass_vel, Particle &ptc, Args &... ptcs) {
  static_assert(calc::all(std::is_same_v<Args, Particle>...), "Type of the 1st argument and the rest should be same!");
  auto dp = centre_mass_pos - COM_p(ptc, ptcs...);
  auto dv = centre_mass_vel - COM_v(ptc, ptcs...);
  ((ptc.pos += dp, ptc.vel += dv), ..., (ptcs.pos += dp, ptcs.vel += dv));
}

template <typename Scalar, typename Particle, typename... Args>
void move_particles(OrbitArgs<Scalar> const &args, Particle &ptc, Args &... ptcs) {
  static_assert(calc::all(std::is_same_v<Args, Particle>...), "Type of the 2nd argument and the rest should be same!");
  auto [cm_pos, cm_vel] = orbit_to_coord(args);
  move_particles(cm_pos, cm_vel, ptc, ptcs...);
}

template <typename Particle, typename... Args>
void move_to_COM_frame(Particle &ptc, Args &... ptcs) {
  using Vector = decltype(COM_v(ptc, ptcs...));
  move_particles(Vector{0, 0, 0}, Vector{0, 0, 0}, ptc, ptcs...);
}

template <typename Particle1, typename Particle2>
inline auto calc_eccentricity(Particle1 const &p1, Particle2 const &p2) {
  auto m1 = M_tot(p1);
  auto m2 = M_tot(p2);
  auto dp = COM_p(p1) - COM_p(p2);
  auto dv = COM_v(p1) - COM_v(p2);
  return calc_eccentricity(consts::G * (m1 + m2), dp, dv);
}

template <typename Particle1, typename Particle2>
inline auto calc_semi_major_axis(Particle1 const &p1, Particle2 const &p2) {
  auto m1 = M_tot(p1);
  auto m2 = M_tot(p2);
  auto dp = COM_p(p1) - COM_p(p2);
  auto dv = COM_v(p1) - COM_v(p2);
  return calc_semi_major_axis(consts::G * (m1 + m2), dp, dv);
}

template <typename Particle1, typename Particle2>
inline auto calc_a_e(Particle1 const &p1, Particle2 const &p2) {
  auto m1 = M_tot(p1);
  auto m2 = M_tot(p2);
  auto dp = COM_p(p1) - COM_p(p2);
  auto dv = COM_v(p1) - COM_v(p2);
  return calc_a_e(consts::G * (m1 + m2), dp, dv);
}

template <typename Particle1, typename Particle2>
inline auto period(Particle1 const &p1, Particle2 const &p2) {
  auto m1 = M_tot(p1);
  auto m2 = M_tot(p2);
  return period(m1, m2, calc_semi_major_axis(p1, p2));
}

template <typename Cluster>
auto E_k(Cluster &&ptc) {
  if constexpr (is_container_v<Cluster>) {
    using Scalar = typename Cluster::value_type::Scalar;
    using Vector = typename Cluster::value_type::Vector;
    Scalar kinetic_energy = 0;

    for (auto &p : ptc) {
      kinetic_energy += p.mass * dot(p.vel, p.vel);
    }

    return 0.5 * kinetic_energy;
  } else {
    return 0.5 * ptc.mass * dot(ptc.vel, ptc.vel);
  }
}

template <typename Particle, typename... Args>
auto E_k(Particle &&ptc, Args &&... args) {
  static_assert(calc::all(std::is_same_v<Args, Particle>...), "Type of the 1st argument and the rest should be same!");

  return 0.5 * ((args.mass * dot(args.vel, args.vel)) + ... + (ptc.mass * dot(ptc.vel, ptc.vel)));
}

template <typename Cluster>
auto E_p(Cluster &&ptc) {
  if constexpr (is_container_v<Cluster>) {
    typename Cluster::value_type::Scaler potential_energy = 0;

    for (auto i = ptc.begin(); i < ptc.end(); ++i) {
      for (auto j = i + 1; j < ptc.end(); ++j) {
        potential_energy -= consts::G * i->mass * j->mass / distance(i->pos, j->pos);
      }
    }
    return potential_energy;
  } else {
    return static_cast<typename Cluster::Scaler>(0);
  }
}

template <typename Particle, typename... Args>
auto E_p(Particle &&ptc, Args &&... args) {
  static_assert(calc::all(std::is_same_v<Args, Particle>...), "Type of the 1st argument and the rest should be same!");
  return E_p(cluster(ptc, args...));
}

CREATE_MEMBER_CHECK(mass);
CREATE_MEMBER_CHECK(radius);

template <typename Particle, typename... Args>
auto E_tot(Particle &&ptc, Args &&... args) {
  return E_k(std::forward<Particle>(ptc), std::forward<Args>(args)...) +
         E_p(std::forward<Particle>(ptc), std::forward<Args>(args)...);
}

template <typename Cluster>
auto cluster_size(Cluster &&ptc) {
  if constexpr (is_container_v<Cluster>) {
    typename Cluster::value_type::Scaler R_max = 0;

    if (ptc.begin() != ptc.end()) {
      R_max = cluster_size(*(ptc.begin()));
    }

    for (auto i = ptc.begin(); i < ptc.end(); ++i) {
      for (auto j = i + 1; j < ptc.end(); ++j) {
        auto r = distance(i->pos, j->pos);
        if (r > R_max) {
          R_max = r;
        }
      }
    }
    return R_max;
  } else {
    if constexpr (HAS_MEMBER(Cluster, radius)) {
      return ptc.radius;
    } else {
      return static_cast<typename Cluster::Scaler>(0);
    }
  }
}

template <typename Particle, typename... Args>
auto cluster_size(Particle &&ptc, Args &&... args) {
  static_assert(calc::all(std::is_same_v<Args, Particle>...), "Type of the 1st argument and the rest should be same!");
  return cluster_size(cluster(ptc, args...));
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
#endif  // SPACEHUB_OPERATIONS_HPP
