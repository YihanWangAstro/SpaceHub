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
 * @file archain.hpp
 *
 * Header file.
 */
#ifndef SPACEHUB_ARCHAIN_HPP
#define SPACEHUB_ARCHAIN_HPP

#include "chain.hpp"
#include "regu-system.hpp"

namespace space::particle_system {

/*---------------------------------------------------------------------------*\
    Class ARchainSystem Declaration
\*---------------------------------------------------------------------------*/
/**
 *
 * @tparam Particles
 * @tparam Interactions
 * @tparam RegType
 */
template <typename Particles, typename Interactions, ReguType RegType>
class ARchainSystem : public ParticleSystem<ARchainSystem<Particles, Interactions, RegType>> {
 public:
  // Type members
  SPACEHUB_USING_TYPE_SYSTEM_OF(Particles);

  using Base = ParticleSystem<ARchainSystem<Particles, Interactions, RegType>>;

  using Particle = typename Particles::Particle;

  // Static public members
  static constexpr ReguType regu_type{RegType};

  // Constructors
  SPACEHUB_MAKE_CONSTRUCTORS(ARchainSystem, delete, default, default, default, default);

  template <typename STL>
  ARchainSystem(Scalar time, STL const &particle_set);

  // Public methods
  SPACEHUB_STD_ACCESSOR(auto, chain_pos, chain_pos_);

  SPACEHUB_STD_ACCESSOR(auto, chain_vel, chain_vel_);

  SPACEHUB_STD_ACCESSOR(auto, index, index_);

  SPACEHUB_STD_ACCESSOR(auto, omega, regu_.omega());

  SPACEHUB_STD_ACCESSOR(auto, bindE, regu_.bindE());

  SPACEHUB_CONDITIONAL_ACCESSOR(HAS_METHOD(Particles, radius), auto, radius, ptcl_.radius());

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

  void impl_post_iter_process();

  template <typename STL>
  void impl_to_linear_container(STL &stl_ranges);

  template <typename STL>
  void impl_load_from_linear_container(STL const &stl_ranges);

  // Friend functions
  template <typename P, typename F, ReguType R>
  friend std::ostream &operator<<(std::ostream &os, ARchainSystem<P, F, R> const &ps);

  template <typename P, typename F, ReguType R>
  friend std::istream &operator>>(std::istream &is, ARchainSystem<P, F, R> &ps);

 private:
  // Private methods
  void chain_advance(Coord &var, Coord &chain_var, Coord const &chain_increment, Scalar phy_time);

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

  Coord chain_pos_;

  Coord chain_vel_;

  Coord chain_acc_;

  IdxArray index_;

  IdxArray new_index_;

  std::conditional_t<Interactions::ext_vel_dep, Coord, Empty> aux_vel_;

  std::conditional_t<Interactions::ext_vel_dep, Coord, Empty> chain_aux_vel_;
};

/*---------------------------------------------------------------------------*\
    Class ARchainSystem Implementation
\*---------------------------------------------------------------------------*/
template <typename Particles, typename Interactions, ReguType RegType>
template <typename STL>
ARchainSystem<Particles, Interactions, RegType>::ARchainSystem(Scalar time, const STL &particle_set)
    : ptcl_(time, particle_set),
      accels_(particle_set.size()),
      regu_(ptcl_),
      chain_pos_(particle_set.size()),
      chain_vel_(particle_set.size()),
      chain_acc_(particle_set.size()),
      index_(particle_set.size()),
      new_index_(particle_set.size()) {
  static_assert(is_ranges_v<STL>, "Only STL-like container can be used");
  Chain::calc_chain_index(ptcl_.pos(), index_);
  Chain::calc_chain(ptcl_.pos(), chain_pos(), index());
  Chain::calc_chain(ptcl_.vel(), chain_vel(), index());
  if constexpr (Interactions::ext_vel_dep) {
    aux_vel_ = ptcl_.vel();
    chain_aux_vel_ = chain_vel_;
  }
}

template <typename Particles, typename Interactions, ReguType RegType>
size_t ARchainSystem<Particles, Interactions, RegType>::impl_number() const {
  return ptcl_.number();
}

template <typename Particles, typename Interactions, ReguType RegType>
void ARchainSystem<Particles, Interactions, RegType>::impl_advance_time(Scalar step_size) {
  Scalar phy_time = regu_.eval_pos_phy_time(*this, step_size);
  ptcl_.time() += phy_time;
}

template <typename Particles, typename Interactions, ReguType RegType>
void ARchainSystem<Particles, Interactions, RegType>::impl_advance_pos(Scalar step_size, const Coord &velocity) {
  Scalar phy_time = regu_.eval_pos_phy_time(*this, step_size);
  Chain::calc_chain(velocity, chain_acc_, index());//borrow chain_acc_ as chain velocity buffer.
  chain_advance(ptcl_.pos(), chain_pos(), chain_acc_, phy_time);
}

template <typename Particles, typename Interactions, ReguType RegType>
void ARchainSystem<Particles, Interactions, RegType>::impl_advance_vel(Scalar step_size, const Coord &acceleration) {
  Scalar phy_time = regu_.eval_vel_phy_time(*this, step_size);
  Chain::calc_chain(acceleration, chain_acc_, index());
  chain_advance(ptcl_.vel(), chain_vel(), chain_acc_, phy_time);
}

template <typename Particles, typename Interactions, ReguType RegType>
void ARchainSystem<Particles, Interactions, RegType>::impl_evaluate_acc(Coord &acceleration) const {
  interactions_.eval_acc(*this, acceleration);
}

template <typename Particles, typename Interactions, ReguType RegType>
void ARchainSystem<Particles, Interactions, RegType>::impl_drift(Scalar step_size) {
  Scalar phy_time = regu_.eval_pos_phy_time(*this, step_size);
  chain_advance(ptcl_.pos(), chain_pos(), chain_vel(), phy_time);
  ptcl_.time() += phy_time;
}

template <typename Particles, typename Interactions, ReguType RegType>
void ARchainSystem<Particles, Interactions, RegType>::impl_kick(Scalar step_size) {
  Scalar phy_time = regu_.eval_vel_phy_time(*this, step_size);
  Scalar half_time = 0.5 * phy_time;
  eval_vel_indep_acc();

  if constexpr (Interactions::ext_vel_dep) {
    kick_pseu_vel(half_time);
    kick_real_vel(phy_time);
    kick_pseu_vel(half_time);
  } else {
    /*Chain::calc_chain(accels_.tot_vel_indep_acc(), chain_acc_, index());
    chain_advance(ptcl_.vel(), chain_vel(), chain_acc_, half_time);
    advance_omega(ptcl_.vel(), accels_.newtonian_acc(), phy_time);
    if constexpr (Interactions::ext_vel_indep) {
      advance_bindE(ptcl_.vel(), accels_.ext_vel_indep_acc(), phy_time);
    }
    chain_advance(ptcl_.vel(), chain_vel(), chain_acc_, half_time);*/

    advance_omega(ptcl_.vel(), accels_.newtonian_acc(), half_time);
    if constexpr (Interactions::ext_vel_indep) {
      advance_bindE(ptcl_.vel(), accels_.ext_vel_indep_acc(), half_time);
    }

    Chain::calc_chain(accels_.tot_vel_indep_acc(), chain_acc_, index());
    chain_advance(ptcl_.vel(), chain_vel(), chain_acc_, phy_time);

    advance_omega(ptcl_.vel(), accels_.newtonian_acc(), half_time);
    if constexpr (Interactions::ext_vel_indep) {
      advance_bindE(ptcl_.vel(), accels_.ext_vel_indep_acc(), half_time);
    }
  }
}

template <typename Particles, typename Interactions, ReguType RegType>
void ARchainSystem<Particles, Interactions, RegType>::impl_pre_iter_process() {
  if constexpr (Interactions::ext_vel_dep) {
    aux_vel_ = ptcl_.vel();
    chain_aux_vel_ = chain_vel_;
  }
}

template <typename Particles, typename Interactions, ReguType RegType>
void ARchainSystem<Particles, Interactions, RegType>::impl_post_iter_process() {
  Chain::calc_chain_index(ptcl_.pos(), new_index_);
  if (new_index_ != index_) {
    Chain::update_chain(chain_pos_, index_, new_index_);
    Chain::calc_cartesian(ptcl_.mass(), chain_pos_, ptcl_.pos(), new_index_);
    Chain::update_chain(chain_vel_, index_, new_index_);
    Chain::calc_cartesian(ptcl_.mass(), chain_vel_, ptcl_.vel(), new_index_);
    index_ = new_index_;
  }
}

template <typename Particles, typename Interactions, ReguType RegType>
template <typename STL>
void ARchainSystem<Particles, Interactions, RegType>::impl_to_linear_container(STL &stl_ranges) {
  stl_ranges.clear();
  stl_ranges.reserve(impl_number() * 6 + 3);
  stl_ranges.emplace_back(impl_time());
  stl_ranges.emplace_back(omega());
  stl_ranges.emplace_back(bindE());
  add_coords_to(stl_ranges, chain_pos_);
  add_coords_to(stl_ranges, chain_vel_);

    /*stl_ranges.reserve(impl_number() * 12 + 3);
    stl_ranges.emplace_back(impl_time());
    stl_ranges.emplace_back(omega());
    stl_ranges.emplace_back(bindE());
    add_coords_to(stl_ranges, chain_pos_);
    add_coords_to(stl_ranges, chain_vel_);
    add_coords_to(stl_ranges, impl_pos());
    add_coords_to(stl_ranges, impl_vel());*/
}

template <typename Particles, typename Interactions, ReguType RegType>
template <typename STL>
void ARchainSystem<Particles, Interactions, RegType>::impl_load_from_linear_container(const STL &stl_ranges) {
  auto begin = stl_ranges.begin();
  impl_time() = *begin;
  omega() = *(begin + 1);
  bindE() = *(begin + 2);

  size_t len = impl_number() * 3;

  auto pos_begin = begin + 3;
  auto pos_end = pos_begin + len;
  auto vel_begin = pos_end;
  auto vel_end = vel_begin + len;
  load_to_coords(pos_begin, pos_end, chain_pos_);
  load_to_coords(vel_begin, vel_end, chain_vel_);

  /*auto c_pos_begin = vel_end;
  auto c_pos_end = c_pos_begin + len;
  auto c_vel_begin = c_pos_end;
  auto c_vel_end = c_vel_begin + len;

  load_to_coords(c_pos_begin, c_pos_end, impl_pos());
  load_to_coords(c_vel_begin, c_vel_end, impl_vel());*/

  Chain::calc_cartesian(ptcl_.mass(), chain_pos_, impl_pos(), index());
  Chain::calc_cartesian(ptcl_.mass(), chain_vel_, impl_vel(), index());
}

template <typename Particles, typename Interactions, ReguType RegType>
std::ostream &operator<<(std::ostream &os, ARchainSystem<Particles, Interactions, RegType> const &ps) {
  os << ps.ptcl_;
  return os;
}

template <typename Particles, typename Interactions, ReguType RegType>
std::istream &operator>>(std::istream &is, ARchainSystem<Particles, Interactions, RegType> &ps) {
  is >> ps.ptcl_;
  return is;
}

template <typename Particles, typename Interactions, ReguType RegType>
void ARchainSystem<Particles, Interactions, RegType>::chain_advance(Coord &var, Coord &chain_var,
                                                                    const Coord &chain_increment, Scalar phy_time) {
  calc::coord_advance(chain_var, chain_increment, phy_time);
  Chain::calc_cartesian(ptcl_.mass(), chain_var, var, index());
}

template <typename Particles, typename Interactions, ReguType RegType>
void ARchainSystem<Particles, Interactions, RegType>::eval_vel_indep_acc() {
  interactions_.eval_newtonian_acc(*this, accels_.newtonian_acc());

  if constexpr (Interactions::ext_vel_indep) {
    interactions_.eval_extra_vel_indep_acc(*this, accels_.ext_vel_indep_acc());
    calc::coord_add(accels_.tot_vel_indep_acc(), accels_.ext_vel_indep_acc(), accels_.newtonian_acc());
  } else {
    accels_.tot_vel_indep_acc() = accels_.newtonian_acc();
  }
}

template <typename Particles, typename Interactions, ReguType RegType>
void ARchainSystem<Particles, Interactions, RegType>::advance_omega(const Coord &velocity, const Coord &d_omega_dr,
                                                                    Scalar phy_time) {
  if constexpr(regu_type == RegType::TTL) {
    Scalar d_omega = calc::coord_contract_to_scalar(ptcl_.mass(), velocity, d_omega_dr);
    regu_.omega() += d_omega * phy_time;
  }
}

template <typename Particles, typename Interactions, ReguType RegType>
void ARchainSystem<Particles, Interactions, RegType>::advance_bindE(const Coord &velocity, const Coord &d_bindE_dr,
                                                                    Scalar phy_time) {
  if constexpr ((Interactions::ext_vel_indep || Interactions::ext_vel_dep) && regu_type == RegType::LogH ) {
    Scalar d_bindE = -calc::coord_contract_to_scalar(ptcl_.mass(), velocity, d_bindE_dr);
    regu_.bindE() += d_bindE * phy_time;
  }
}

template <typename Particles, typename Interactions, ReguType RegType>
void ARchainSystem<Particles, Interactions, RegType>::kick_pseu_vel(Scalar phy_time) {
  interactions_.eval_extra_vel_dep_acc(*this, accels_.ext_vel_dep_acc());
  calc::coord_add(accels_.acc(), accels_.tot_vel_indep_acc(), accels_.ext_vel_dep_acc());
  Chain::calc_chain(accels_.acc(), chain_acc_, index());
  chain_advance(aux_vel_, chain_aux_vel_, chain_acc_, phy_time);
}

template <typename Particles, typename Interactions, ReguType RegType>
void ARchainSystem<Particles, Interactions, RegType>::kick_real_vel(Scalar phy_time) {
  std::swap(aux_vel_, ptcl_.vel());
  std::swap(chain_aux_vel_, chain_vel());
  interactions_.eval_extra_vel_dep_acc(*this, accels_.ext_vel_dep_acc());
  std::swap(aux_vel_, ptcl_.vel());
  std::swap(chain_aux_vel_, chain_vel());
  calc::coord_add(accels_.acc(), accels_.tot_vel_indep_acc(), accels_.ext_vel_dep_acc());

  Chain::calc_chain(accels_.acc(), chain_acc_, index());
  chain_advance(ptcl_.vel(), chain_vel(), chain_acc_, phy_time);

  advance_omega(aux_vel_, accels_.newtonian_acc(), phy_time);

  if constexpr (Interactions::ext_vel_indep) {
    calc::coord_add(accels_.acc(), accels_.ext_vel_indep_acc(), accels_.ext_vel_dep_acc());
    advance_bindE(aux_vel_, accels_.acc(), phy_time);
  } else {
    advance_bindE(aux_vel_, accels_.ext_vel_dep_acc(), phy_time);
  }
}
}  // namespace space::particle_system
#endif  // SPACEHUB_ARCHAIN_HPP
