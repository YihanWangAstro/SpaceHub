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
 * @file regu-system.hpp
 *
 * Header file.
 */
#ifndef SPACEHUB_REGU_SYSTEM_HPP
#define SPACEHUB_REGU_SYSTEM_HPP

#include <type_traits>
#include "../core-computation.hpp"
#include "../interaction/accelerations.hpp"
#include "particle-system.hpp"

namespace space::particle_system {
/**
 *
 */
enum class ReguType { LogH, TTL, None };

/*---------------------------------------------------------------------------*\
    Class Regularization Declaration
\*---------------------------------------------------------------------------*/
/**
 *
 * @tparam Scalar
 * @tparam Type
 */
template <typename Scalar, ReguType Type = ReguType::LogH>
class Regularization {
 public:
  // Constructors
  template <typename Particles>
  explicit Regularization(Particles const &particles);

  // Public methods
  SPACEHUB_STD_ACCESSOR(auto, omega, omega_);

  SPACEHUB_STD_ACCESSOR(auto, bindE, bindE_);

  template <typename Particles>
  Scalar eval_pos_phy_time(Particles const &particles, Scalar step_size) const;

  template <typename Particles>
  Scalar eval_vel_phy_time(Particles const &particles, Scalar step_size) const;

 private:
  // Private methods
  template <typename Particles>
  inline Scalar capital_omega(Particles const &particles) const;

  // Private members
  Scalar omega_;

  Scalar bindE_;
};

/*---------------------------------------------------------------------------*\
    Class RegularizedSystem Declaration
\*---------------------------------------------------------------------------*/

/**
 * @brief Regularized particle System.
 *
 * Regularied particle system.
 * See details [Explicit Symplectic Algorithms For Time‐Transformed
 * Hamiltonians](https://link.springer.com/article/10.1023%2FA%3A1008368322547) , [A Class of Symplectic Integrators
 * with Adaptive Time Step for Separable Hamiltonian Systems](http://iopscience.iop.org/article/10.1086/301102/meta) and
 * [A Time-Transformed Leapfrog Scheme](https://link.springer.com/article/10.1023%2FA%3A1021149313347).
 * @tparam Particles
 * @tparam Interactions
 */
template <typename Particles, typename Interactions, ReguType RegType>
class RegularizedSystem : public ParticleSystem<RegularizedSystem<Particles, Interactions, RegType>> {
 public:
  // Type members
  SPACEHUB_USING_TYPE_SYSTEM_OF(Particles);

  using Base = ParticleSystem<RegularizedSystem<Particles, Interactions, RegType>>;

  using Particle = typename Particles::Particle;

  // Constructors
  SPACEHUB_MAKE_CONSTRUCTORS(RegularizedSystem, delete, default, default, default, default);

  template <typename STL>
  RegularizedSystem(Scalar time, STL const &particle_set);

  // Static members
  static constexpr ReguType regu_type{RegType};

  // Public methods
  SPACEHUB_STD_ACCESSOR(auto, omega, regu_.omega());

  SPACEHUB_STD_ACCESSOR(auto, bindE, regu_.bindE());

  // Friend functions
  template <typename P, typename F, ReguType R>
  friend std::ostream &operator<<(std::ostream &os, RegularizedSystem<P, F, R> const &ps);

  template <typename P, typename F, ReguType R>
  friend std::istream &operator>>(std::istream &is, RegularizedSystem<P, F, R> &ps);

  CRTP_IMPL :
      // CRTP implementation
      SPACEHUB_STD_ACCESSOR(auto, impl_mass, ptcl_.mass());

  SPACEHUB_STD_ACCESSOR(auto, impl_idn, ptcl_.idn());

  SPACEHUB_STD_ACCESSOR(auto, impl_pos, ptcl_.pos());

  SPACEHUB_STD_ACCESSOR(auto, impl_vel, ptcl_.vel());

  SPACEHUB_STD_ACCESSOR(auto, impl_time, ptcl_.time());

  [[nodiscard]] size_t impl_number() const;

  void impl_advance_time(Scalar step_size);

  void impl_advance_pos(Scalar step_size, Coord const &velocity);

  void impl_advance_vel(Scalar step_size, Coord const &acceleration);

  void impl_evaluate_acc(Coord &acceleration) const;

  void impl_drift(Scalar step_size);

  void impl_kick(Scalar step_size);

  void impl_pre_iter_process();

  template <typename STL>
  void impl_to_linear_container(STL &stl_ranges);

  template <typename STL>
  void impl_load_from_linear_container(STL const &stl_ranges);

 private:
  // Private methods
  void eval_vel_indep_acc();

  void advance_omega(Coord const &velocity, Coord const &d_omega_dr, Scalar phy_time);

  void advance_bindE(Coord const &velocity, Coord const &d_bindE_dr, Scalar phy_time);

  void kick_pseu_vel(Scalar phy_time);

  void kick_real_vel(Scalar phy_time);

  // Private members
  Particles ptcl_;
  Interactions interactions_;
  interactions::Accelerations<Interactions, Coord> accels_;
  Regularization<Scalar, RegType> regu_;
  std::conditional_t<Interactions::ext_vel_dep, Coord, Empty> aux_vel_;
};

/*---------------------------------------------------------------------------*\
    Class RegularizedSystem Implementation
\*---------------------------------------------------------------------------*/
template <typename Particles, typename Interactions, ReguType RegType>
template <typename STL>
RegularizedSystem<Particles, Interactions, RegType>::RegularizedSystem(Scalar time, const STL &particle_set)
    : ptcl_(time, particle_set), accels_(particle_set.size()), regu_(ptcl_) {
  static_assert(is_ranges_v<STL>, "Only STL-like container can be used");
  if constexpr (Interactions::ext_vel_dep) {
    aux_vel_ = ptcl_.vel();
  }
}

template <typename Particles, typename Interactions, ReguType RegType>
std::istream &operator>>(std::istream &is, RegularizedSystem<Particles, Interactions, RegType> &ps) {
  is >> ps.ptcl_;
  return is;
}

template <typename Particles, typename Interactions, ReguType RegType>
std::ostream &operator<<(std::ostream &os, const RegularizedSystem<Particles, Interactions, RegType> &ps) {
  os << ps.ptcl_;
  return os;
}

template <typename Particles, typename Interactions, ReguType RegType>
size_t RegularizedSystem<Particles, Interactions, RegType>::impl_number() const {
  return ptcl_.number();
}

template <typename Particles, typename Interactions, ReguType RegType>
void RegularizedSystem<Particles, Interactions, RegType>::impl_advance_time(Scalar step_size) {
  Scalar phy_time = regu_.eval_pos_phy_time(ptcl_, step_size);
  ptcl_.time() += phy_time;
}

template <typename Particles, typename Interactions, ReguType RegType>
void RegularizedSystem<Particles, Interactions, RegType>::impl_advance_pos(Scalar step_size, const Coord &velocity) {
  Scalar phy_time = regu_.eval_pos_phy_time(ptcl_, step_size);
  calc::coord_advance(ptcl_.pos(), velocity, phy_time);
}

template <typename Particles, typename Interactions, ReguType RegType>
void RegularizedSystem<Particles, Interactions, RegType>::impl_advance_vel(Scalar step_size,
                                                                           const Coord &acceleration) {
  Scalar phy_time = regu_.eval_vel_phy_time(ptcl_, step_size);
  Scalar half_time = 0.5 * phy_time;
  calc::coord_advance(ptcl_.vel(), acceleration, half_time);
  advance_omega(ptcl_.vel(), acceleration, phy_time);
  // TODO : evolve bindE
  calc::coord_advance(ptcl_.vel(), acceleration, half_time);
}

template <typename Particles, typename Interactions, ReguType RegType>
void RegularizedSystem<Particles, Interactions, RegType>::impl_evaluate_acc(Coord &acceleration) const {
  interactions_.eval_acc(ptcl_, acceleration);
}

template <typename Particles, typename Interactions, ReguType RegType>
void RegularizedSystem<Particles, Interactions, RegType>::impl_drift(Scalar step_size) {
  Scalar phy_time = regu_.eval_pos_phy_time(ptcl_, step_size);
  calc::coord_advance(ptcl_.pos(), ptcl_.vel(), phy_time);
  ptcl_.time() += phy_time;
}

template <typename Particles, typename Interactions, ReguType RegType>
void RegularizedSystem<Particles, Interactions, RegType>::impl_kick(Scalar step_size) {
  Scalar phy_time = regu_.eval_vel_phy_time(ptcl_, step_size);
  Scalar half_time = 0.5 * phy_time;

  eval_vel_indep_acc();

  if constexpr (Interactions::ext_vel_dep) {
    kick_pseu_vel(half_time);
    kick_real_vel(phy_time);
    kick_pseu_vel(half_time);
  } else {
    calc::coord_advance(ptcl_.vel(), accels_.tot_vel_indep_acc(), half_time);
    advance_omega(ptcl_.vel(), accels_.newtonian_acc(), phy_time);
    if constexpr (Interactions::ext_vel_indep) {
      advance_bindE(ptcl_.vel(), accels_.ext_vel_indep_acc(), phy_time);
    }
    calc::coord_advance(ptcl_.vel(), accels_.tot_vel_indep_acc(), half_time);
    /*advance_omega(ptcl_.vel(), accels_.newtonian_acc(), half_time);
    if constexpr (Interactions::ext_vel_indep) {
      advance_bindE(ptcl_.vel(), accels_.ext_vel_indep_acc(), half_time);
    }
    calc::coord_advance(ptcl_.vel(), accels_.tot_vel_indep_acc(), phy_time);
    advance_omega(ptcl_.vel(), accels_.newtonian_acc(), half_time);
    if constexpr (Interactions::ext_vel_indep) {
      advance_bindE(ptcl_.vel(), accels_.ext_vel_indep_acc(), half_time);
    }*/
  }
}

template <typename Particles, typename Interactions, ReguType RegType>
void RegularizedSystem<Particles, Interactions, RegType>::impl_pre_iter_process() {
  if constexpr (Interactions::ext_vel_dep) {
    aux_vel_ = ptcl_.vel();
  }
}

template <typename Particles, typename Interactions, ReguType RegType>
template <typename STL>
void RegularizedSystem<Particles, Interactions, RegType>::impl_to_linear_container(STL &stl_ranges) {
  stl_ranges.clear();
  stl_ranges.reserve(impl_number() * 6 + 3);
  stl_ranges.emplace_back(impl_time());
  stl_ranges.emplace_back(omega());
  stl_ranges.emplace_back(bindE());
  add_coords_to(stl_ranges, impl_pos());
  add_coords_to(stl_ranges, impl_vel());
}

template <typename Particles, typename Interactions, ReguType RegType>
template <typename STL>
void RegularizedSystem<Particles, Interactions, RegType>::impl_load_from_linear_container(const STL &stl_ranges) {
  auto begin = stl_ranges.begin();
  impl_time() = *begin;
  omega() = *(begin + 1);
  bindE() = *(begin + 2);

  size_t len = impl_number() * 3;

  auto pos_begin = begin + 3;
  auto pos_end = pos_begin + len;
  auto vel_begin = pos_end;
  auto vel_end = vel_begin + len;
  load_to_coords(pos_begin, pos_end, impl_pos());
  load_to_coords(vel_begin, vel_end, impl_vel());
}

template <typename Particles, typename Interactions, ReguType RegType>
void RegularizedSystem<Particles, Interactions, RegType>::eval_vel_indep_acc() {
  interactions_.eval_newtonian_acc(ptcl_, accels_.newtonian_acc());
  accels_.tot_vel_indep_acc() = accels_.newtonian_acc();
  if constexpr (Interactions::ext_vel_indep) {
    interactions_.eval_extra_vel_indep_acc(ptcl_, accels_.ext_vel_indep_acc());
    calc::coord_add(accels_.tot_vel_indep_acc(), accels_.ext_vel_indep_acc(), accels_.newtonian_acc());
  }
}

template <typename Particles, typename Interactions, ReguType RegType>
void RegularizedSystem<Particles, Interactions, RegType>::advance_omega(const Coord &velocity, const Coord &d_omega_dr,
                                                                        Scalar phy_time) {
  Scalar d_omega = calc::coord_contract_to_scalar(ptcl_.mass(), velocity, d_omega_dr);
  regu_.omega() += d_omega * phy_time;
}

template <typename Particles, typename Interactions, ReguType RegType>
void RegularizedSystem<Particles, Interactions, RegType>::advance_bindE(const Coord &velocity, const Coord &d_bindE_dr,
                                                                        Scalar phy_time) {
  Scalar d_bindE = -calc::coord_contract_to_scalar(ptcl_.mass(), velocity, d_bindE_dr);
  regu_.bindE() += d_bindE * phy_time;
}

template <typename Particles, typename Interactions, ReguType RegType>
void RegularizedSystem<Particles, Interactions, RegType>::kick_pseu_vel(Scalar phy_time) {
  interactions_.eval_extra_vel_dep_acc(ptcl_, accels_.ext_vel_dep_acc());
  calc::coord_add(accels_.acc(), accels_.tot_vel_indep_acc(), accels_.ext_vel_dep_acc());
  calc::coord_advance(aux_vel_, accels_.acc(), phy_time);
}

template <typename Particles, typename Interactions, ReguType RegType>
void RegularizedSystem<Particles, Interactions, RegType>::kick_real_vel(Scalar phy_time) {
  std::swap(aux_vel_, ptcl_.vel());
  interactions_.eval_extra_vel_dep_acc(ptcl_, accels_.ext_vel_dep_acc());
  std::swap(aux_vel_, ptcl_.vel());

  calc::coord_add(accels_.acc(), accels_.tot_vel_indep_acc(), accels_.ext_vel_dep_acc());
  calc::coord_advance(ptcl_.vel(), accels_.acc(), phy_time);

  advance_omega(aux_vel_, accels_.newtonian_acc(), phy_time);

  if constexpr (Interactions::ext_vel_indep) {
    calc::coord_add(accels_.acc(), accels_.ext_vel_indep_acc(), accels_.ext_vel_dep_acc());
    advance_bindE(aux_vel_, accels_.acc(), phy_time);
  } else {
    advance_bindE(aux_vel_, accels_.ext_vel_dep_acc(), phy_time);
  }
}

/*---------------------------------------------------------------------------*\
    Class Regularization Implementation
\*---------------------------------------------------------------------------*/
template <typename Scalar, ReguType Type>
template <typename Particles>
Regularization<Scalar, Type>::Regularization(const Particles &particles) {
  omega_ = capital_omega(particles);
  bindE_ = -calc::calc_total_energy(particles);
}

template <typename Scalar, ReguType Type>
template <typename Particles>
Scalar Regularization<Scalar, Type>::eval_pos_phy_time(const Particles &particles, Scalar step_size) const {
  if constexpr (Type == ReguType::LogH) {
    return step_size / (bindE_ + calc::calc_kinetic_energy(particles));
  } else if constexpr (Type == ReguType::TTL) {
    return step_size / omega_;
  } else if constexpr (Type == ReguType::None) {
    return step_size;
  } else {
    spacehub_abort("Undefined regularization type!");
  }
}

template <typename Scalar, ReguType Type>
template <typename Particles>
Scalar Regularization<Scalar, Type>::eval_vel_phy_time(const Particles &particles, Scalar step_size) const {
  if constexpr (Type == ReguType::LogH) {
    return step_size / -calc::calc_potential_energy(particles);
  } else if constexpr (Type == ReguType::TTL) {
    return step_size / capital_omega(particles);
  } else if constexpr (Type == ReguType::None) {
    return step_size;
  } else {
    spacehub_abort("Undefined regularization type!");
  }
}

template <typename Scalar, ReguType Type>
template <typename Particles>
Scalar Regularization<Scalar, Type>::capital_omega(const Particles &particles) const {
  return -calc::calc_potential_energy(particles);
}
}  // namespace space::particle_system

#endif
