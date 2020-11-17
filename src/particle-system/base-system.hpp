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
 * @file base-system.hpp
 *
 * Header file.
 */
#pragma once

#include <type_traits>

#include "../core-computation.hpp"
#include "../interaction/accelerations.hpp"
#include "../spacehub-concepts.hpp"

namespace space::particle_system {

/*---------------------------------------------------------------------------*\
    Class SimpleSystem Declaration
\*---------------------------------------------------------------------------*/
/**
 * @tparam Particles
 * @tparam Interactions
 */
template <concepts::Particles Particles, concepts::Interaction Interactions>
class SimpleSystem {
 public:
  // Type members
  SPACEHUB_USING_TYPE_SYSTEM_OF(Particles);

  using Particle = typename Particles::Particle;

  // Constructors
  SPACEHUB_MAKE_CONSTRUCTORS(SimpleSystem, delete, default, default, default, default);

  /**
   *
   * @tparam STL
   * @param time
   * @param particle_set
   */
  template <typename STL>
  SimpleSystem(Scalar time, STL const &particle_set);

  // Public methods
  SPACEHUB_READ_ACCESSOR(Particles, particles, ptcl_);

  /**
   *
   * @return
   */
  [[nodiscard]] size_t number() const;

  /**
   *
   * @param dt
   */
  void advance_time(Scalar dt);

  /**
   *
   * @param velocity
   * @param step_size
   */
  void advance_pos(Scalar step_size, VectorArray const &velocity);

  /**
   *
   * @param acceleration
   * @param step_size
   */
  void advance_vel(Scalar step_size, VectorArray const &acceleration);

  /**
   *
   * @param acceleration
   */
  void evaluate_acc(VectorArray &acceleration) const;

  /**
   *
   * @param step_size
   */
  void drift(Scalar step_size);

  /**
   *
   * @param step_size
   */
  void kick(Scalar step_size);

  /**
   *
   */
  void pre_iter_process();

  void post_iter_process(){};

  /**
   *
   * @tparam STL
   * @param stl_ranges
   */
  template <typename STL>
  void to_linear_container(STL &stl_ranges);

  /**
   *
   * @tparam STL
   * @param stl_ranges
   */
  template <typename STL>
  void load_from_linear_container(STL const &stl_ranges);

  // Friend functions
  template <typename P, typename F>
  friend std::ostream &operator<<(std::ostream &os, SimpleSystem<P, F> const &ps);

  template <typename P, typename F>
  friend std::istream &operator>>(std::istream &is, SimpleSystem<P, F> &ps);

 private:
  // Private methods
  /**
   *
   */
  void eval_vel_indep_acc();

  /**
   *
   * @param step_size
   */
  void kick_pseu_vel(Scalar step_size);

  /**
   *
   * @param step_size
   */
  void kick_real_vel(Scalar step_size);

  // Private members
  Particles ptcl_;

  Interactions interactions_;

  interactions::Accelerations<Interactions, VectorArray> accels_;

  std::conditional_t<Interactions::ext_vel_dep, VectorArray, Empty> aux_vel_;
};
}  // namespace space::particle_system

namespace space::particle_system {
/*---------------------------------------------------------------------------*\
    Class SimpleSystem Implementation
\*---------------------------------------------------------------------------*/
template <Concept_Particles Particles, Concept_Interactions Interactions>
template <typename STL>
SimpleSystem<Particles, Interactions>::SimpleSystem(Scalar time, const STL &particle_set)
    : ptcl_(time, particle_set), accels_(particle_set.size()) {
  static_assert(is_ranges_v<STL>, "Only STL-like container can be used");
  if constexpr (Interactions::ext_vel_dep) {
    aux_vel_ = ptcl_.vel();
  }
}

template <Concept_Particles Particles, Concept_Interactions Interactions>
template <typename STL>
void SimpleSystem<Particles, Interactions>::load_from_linear_container(const STL &stl_ranges) {
  auto begin = stl_ranges.begin();
  impl_time() = *begin;
  size_t len = impl_number() * 3;
  auto pos_begin = begin + 1;
  auto pos_end = pos_begin + len;
  auto vel_begin = pos_end;
  auto vel_end = vel_begin + len;
  load_to_coords(pos_begin, pos_end, impl_pos());
  load_to_coords(vel_begin, vel_end, impl_vel());
}

template <Concept_Particles Particles, Concept_Interactions Interactions>
template <typename STL>
void SimpleSystem<Particles, Interactions>::to_linear_container(STL &stl_ranges) {
  stl_ranges.clear();
  stl_ranges.reserve(impl_number() * 6 + 1);
  stl_ranges.emplace_back(impl_time());
  add_coords_to(stl_ranges, impl_pos());
  add_coords_to(stl_ranges, impl_vel());
}

template <Concept_Particles Particles, Concept_Interactions Interactions>
void SimpleSystem<Particles, Interactions>::pre_iter_process() {
  if constexpr (Interactions::ext_vel_dep) {
    aux_vel_ = ptcl_.vel();
  }
}

template <Concept_Particles Particles, Concept_Interactions Interactions>
void SimpleSystem<Particles, Interactions>::kick(Scalar step_size) {
  if constexpr (Interactions::ext_vel_dep) {
    Scalar half_step = 0.5 * step_size;
    eval_vel_indep_acc();
    kick_pseu_vel(half_step);
    kick_real_vel(step_size);
    kick_pseu_vel(half_step);
  } else {
    interactions_.eval_acc(ptcl_, accels_.acc());
    calc::coord_advance(ptcl_.vel(), accels_.acc(), step_size);
  }
}

template <Concept_Particles Particles, Concept_Interactions Interactions>
void SimpleSystem<Particles, Interactions>::drift(Scalar step_size) {
  ptcl_.time() += step_size;
  calc::coord_advance(ptcl_.pos(), ptcl_.vel(), step_size);
}

template <Concept_Particles Particles, Concept_Interactions Interactions>
void SimpleSystem<Particles, Interactions>::evaluate_acc(VectorArray &acceleration) const {
  interactions_.eval_acc(ptcl_, acceleration);
}

template <Concept_Particles Particles, Concept_Interactions Interactions>
void SimpleSystem<Particles, Interactions>::advance_vel(Scalar step_size, const VectorArray &acceleration) {
  calc::coord_advance(ptcl_.vel(), acceleration, step_size);
}

template <Concept_Particles Particles, Concept_Interactions Interactions>
void SimpleSystem<Particles, Interactions>::advance_pos(Scalar step_size, const VectorArray &velocity) {
  calc::coord_advance(ptcl_.pos(), velocity, step_size);
}

template <Concept_Particles Particles, Concept_Interactions Interactions>
void SimpleSystem<Particles, Interactions>::advance_time(Scalar dt) {
  ptcl_.time() += dt;
}

template <Concept_Particles Particles, Concept_Interactions Interactions>
size_t SimpleSystem<Particles, Interactions>::number() const {
  return ptcl_.number();
}

template <Concept_Particles Particles, Concept_Interactions Interactions>
void SimpleSystem<Particles, Interactions>::kick_real_vel(Scalar step_size) {
  std::swap(aux_vel_, ptcl_.vel());
  interactions_.eval_extra_vel_dep_acc(ptcl_, accels_.ext_vel_dep_acc());
  std::swap(aux_vel_, ptcl_.vel());
  calc::coord_add(accels_.acc(), accels_.tot_vel_indep_acc(), accels_.ext_vel_dep_acc());
  calc::coord_advance(ptcl_.vel(), accels_.acc(), step_size);
}

template <Concept_Particles Particles, Concept_Interactions Interactions>
void SimpleSystem<Particles, Interactions>::kick_pseu_vel(Scalar step_size) {
  interactions_.eval_extra_vel_dep_acc(ptcl_, accels_.ext_vel_dep_acc());
  calc::coord_add(accels_.acc(), accels_.tot_vel_indep_acc(), accels_.ext_vel_dep_acc());
  calc::coord_advance(aux_vel_, accels_.acc(), step_size);
}

template <Concept_Particles Particles, Concept_Interactions Interactions>
void SimpleSystem<Particles, Interactions>::eval_vel_indep_acc() {
  interactions_.eval_newtonian_acc(ptcl_, accels_.tot_vel_indep_acc());
  if constexpr (Interactions::ext_vel_indep) {
    interactions_.eval_extra_vel_indep_acc(ptcl_, accels_.ext_vel_indep_acc());
    calc::coord_add(accels_.tot_vel_indep_acc(), accels_.tot_vel_indep_acc(), accels_.ext_vel_indep_acc());
  }
}

template <Concept_Particles Particles, Concept_Interactions Interactions>
std::ostream &operator<<(std::ostream &os, SimpleSystem<Particles, Interactions> const &ps) {
  os << ps.ptcl_;
  return os;
}

template <Concept_Particles Particles, Concept_Interactions Interactions>
std::istream &operator>>(std::istream &is, SimpleSystem<Particles, Interactions> &ps) {
  is >> ps.ptcl_;
  return is;
}
}  // namespace space::particle_system
