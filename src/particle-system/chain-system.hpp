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
 * @file chain-system.hpp
 *
 * Header file.
 */
#ifndef SPACEHUB_CHAIN_SYSTEM_HPP
#define SPACEHUB_CHAIN_SYSTEM_HPP

#include "../core-computation.hpp"
#include "particle-system.hpp"
#include "chain.hpp"
#include <type_traits>

namespace space::particle_system {

  /*---------------------------------------------------------------------------*\
      Class ChainSystem Declaration
  \*---------------------------------------------------------------------------*/
  /**
   *
   * @tparam Particles
   * @tparam Interactions
   */
  template<typename Particles, typename Interactions>
  class ChainSystem : public ParticleSystem<ChainSystem<Particles, Interactions>> {
  public:
    //Type members
    SPACEHUB_USING_TYPE_SYSTEM_OF(Particles);

    using Base = ParticleSystem<ChainSystem<Particles, Interactions>>;

    using Particle = typename Particles::Particle;

    //Constructors
    SPACEHUB_MAKE_CONSTRUCTORS(ChainSystem, delete, default, default, default, default);

    template<typename STL>
    ChainSystem(Scalar time, STL const &particle_set);

    //Public methods
    SPACEHUB_STD_ACCESSOR(auto, chain_pos, chain_pos_);

    SPACEHUB_STD_ACCESSOR(auto, chain_vel, chain_vel_);

    SPACEHUB_STD_ACCESSOR(auto, index, index_);

    //Friend functions
    template<typename P, typename F>
    friend std::ostream &operator<<(std::ostream &os, ChainSystem<P, F> const &ps);

    template<typename P, typename F>
    friend std::istream &operator>>(std::istream &is, ChainSystem<P, F> &ps);

    CRTP_IMPL:
    //CRTP implementation
    SPACEHUB_STD_ACCESSOR(auto, impl_mass, ptcl_.mass());

    SPACEHUB_STD_ACCESSOR(auto, impl_idn, ptcl_.idn());

    SPACEHUB_STD_ACCESSOR(auto, impl_pos, ptcl_.pos());

    SPACEHUB_STD_ACCESSOR(auto, impl_vel, ptcl_.vel());

    SPACEHUB_STD_ACCESSOR(auto, impl_time, ptcl_.time());

    [[nodiscard]] size_t impl_number() const;

    void impl_advance_time(Scalar dt);

    void impl_advance_pos(Coord const &velocity, Scalar step_size);

    void impl_advance_vel(Coord const &acceleration, Scalar step_size);

    void impl_evaluate_acc(Coord &acceleration) const;

    void impl_drift(Scalar step_size);

    void impl_kick(Scalar step_size);

    void impl_pre_iter_process();

    void impl_post_iter_process();

    template<typename STL>
    void impl_to_linear_container(STL &stl_ranges);

    template<typename STL>
    void impl_load_from_linear_container(STL const &stl_ranges);

  private:
    //Private methods
    void chain_advance(Coord &var, Coord &chain_var, Coord &chain_increment, Scalar step_size);

    void eval_vel_indep_acc();

    void kick_pseu_vel(Scalar step_size);

    void kick_real_vel(Scalar step_size);

    //Private members
    Particles ptcl_;
    Interactions interactions_;
    interactions::Accelerations <Interactions, Coord> accels_{};
    Coord chain_pos_;
    Coord chain_vel_;
    Coord chain_acc_;

    IdxArray index_;
    IdxArray new_index_;

    std::conditional_t<Interactions::ext_vel_dep, Coord, Empty> aux_vel_;
    std::conditional_t<Interactions::ext_vel_dep, Coord, Empty> chain_aux_vel_;
  };

  /*---------------------------------------------------------------------------*\
      Class ChainSystem Implementation
  \*---------------------------------------------------------------------------*/
  template<typename Particles, typename Interactions>
  template<typename STL>
  ChainSystem<Particles, Interactions>::ChainSystem(Scalar time, const STL &particle_set)
          : ptcl_(time, particle_set),
            chain_pos_(particle_set.size()),
            chain_vel_(particle_set.size()),
            index_(particle_set.size()),
            new_index_(particle_set.size()),
            accels_(particle_set.size()),
            chain_acc_(particle_set.size()) {
    static_assert(is_container_v<STL>, "Only STL-like container can be used");
    Chain::calc_chain_index(ptcl_.pos(), index_);
    Chain::calc_chain(ptcl_.pos(), chain_pos(), index());
    Chain::calc_chain(ptcl_.vel(), chain_vel(), index());

    if constexpr (Interactions::ext_vel_dep) {
      aux_vel_ = ptcl_.vel();
      chain_aux_vel_ = chain_vel_;
    }
  }

  template<typename Particles, typename Interactions>
  size_t ChainSystem<Particles, Interactions>::impl_number() const {
    return ptcl_.number();
  }

  template<typename Particles, typename Interactions>
  void ChainSystem<Particles, Interactions>::impl_advance_time(Scalar dt) {
    ptcl_.time() += dt;
  }

  template<typename Particles, typename Interactions>
  void ChainSystem<Particles, Interactions>::impl_advance_pos(const Coord &velocity, Scalar step_size) {
    Chain::calc_chain(velocity, chain_vel(), index());
    chain_advance(ptcl_.pos(), chain_pos(), chain_vel(), step_size);
  }

  template<typename Particles, typename Interactions>
  void ChainSystem<Particles, Interactions>::impl_advance_vel(const Coord &acceleration, Scalar step_size) {
    Chain::calc_chain(acceleration, chain_acc_, index());
    chain_advance(ptcl_.vel(), chain_vel(), chain_acc_, step_size);
  }

  template<typename Particles, typename Interactions>
  void ChainSystem<Particles, Interactions>::impl_evaluate_acc(Coord &acceleration) const {
    interactions_.eval_acc(*this, acceleration);
  }

  template<typename Particles, typename Interactions>
  void ChainSystem<Particles, Interactions>::impl_drift(Scalar step_size) {
    ptcl_.time() += step_size;
    chain_advance(ptcl_.pos(), chain_pos(), chain_vel(), step_size);
  }

  template<typename Particles, typename Interactions>
  void ChainSystem<Particles, Interactions>::impl_kick(Scalar step_size) {
    if constexpr (Interactions::ext_vel_dep) {
      Scalar half_step = 0.5 * step_size;
      eval_vel_indep_acc();
      kick_pseu_vel(half_step);
      kick_real_vel(step_size);
      kick_pseu_vel(half_step);
    } else {
      interactions_.eval_acc(*this, accels_.acc());
      impl_advance_vel(accels_.acc(), step_size);
    }
  }

  template<typename Particles, typename Interactions>
  void ChainSystem<Particles, Interactions>::impl_pre_iter_process() {
    if constexpr (Interactions::ext_vel_dep) {
      aux_vel_ = ptcl_.vel();
      chain_aux_vel_ = chain_vel_;
    }
  }

  template<typename Particles, typename Interactions>
  void ChainSystem<Particles, Interactions>::impl_post_iter_process() {
    Chain::calc_chain_index(ptcl_.pos(), new_index_);
    if (new_index_ != index_) {
      Chain::update_chain(chain_pos_, index_, new_index_);
      Chain::calc_cartesian(ptcl_.mass(), chain_pos_, ptcl_.pos(), new_index_);
      Chain::update_chain(chain_vel_, index_, new_index_);
      Chain::calc_cartesian(ptcl_.mass(), chain_vel_, ptcl_.vel(), new_index_);
      index_ = new_index_;
    }
  }

  template<typename Particles, typename Interactions>
  template<typename STL>
  void ChainSystem<Particles, Interactions>::impl_to_linear_container(STL &stl_ranges) {
    stl_ranges.clear();
    stl_ranges.reserve(impl_number() * 6 + 1);
    stl_ranges.emplace_back(impl_time());
    add_coords_to(stl_ranges, chain_pos_);
    add_coords_to(stl_ranges, chain_vel_);
  }

  template<typename Particles, typename Interactions>
  template<typename STL>
  void ChainSystem<Particles, Interactions>::impl_load_from_linear_container(const STL &stl_ranges) {
    auto begin = stl_ranges.begin();
    impl_time() = *begin;
    size_t len = impl_number() * 3;
    auto pos_begin = begin + 1;
    auto pos_end = pos_begin + len;
    auto vel_begin = pos_end;
    auto vel_end = vel_begin + len;

    load_to_coords(pos_begin, pos_end, chain_pos_);
    load_to_coords(vel_begin, vel_end, chain_vel_);

    Chain::calc_cartesian(ptcl_.mass(), chain_pos_, ptcl_.pos(), index_);
    Chain::calc_cartesian(ptcl_.mass(), chain_vel_, ptcl_.vel(), index_);
  }

  template<typename Particles, typename Interactions>
  std::istream &operator>>(std::istream &is, ChainSystem<Particles, Interactions> &ps) {
    is >> ps.ptcl_;
    return is;
  }

  template<typename Particles, typename Interactions>
  std::ostream &operator<<(std::ostream &os, const ChainSystem<Particles, Interactions> &ps) {
    os << ps.ptcl_;
    return os;
  }

  template<typename Particles, typename Interactions>
  void ChainSystem<Particles, Interactions>::chain_advance(Coord &var, Coord &chain_var, Coord &chain_increment,
                                                           Scalar step_size) {
    calc::coord_advance(chain_var, chain_increment, step_size);
    Chain::calc_cartesian(ptcl_.mass(), chain_var, var, index());
  }

  template<typename Particles, typename Interactions>
  void ChainSystem<Particles, Interactions>::eval_vel_indep_acc() {
    interactions_.eval_newtonian_acc(*this, accels_.tot_vel_indep_acc());
    if constexpr (Interactions::ext_vel_indep) {
      interactions_.eval_extra_vel_indep_acc(*this, accels_.ext_vel_indep_acc());
      calc::coord_add(accels_.tot_vel_indep_acc(), accels_.tot_vel_indep_acc(), accels_.ext_vel_indep_acc());
    }
  }

  template<typename Particles, typename Interactions>
  void ChainSystem<Particles, Interactions>::kick_pseu_vel(Scalar step_size) {
    interactions_.eval_extra_vel_dep_acc(*this, accels_.ext_vel_dep_acc());
    calc::coord_add(accels_.acc(), accels_.tot_vel_indep_acc(), accels_.ext_vel_dep_acc());
    Chain::calc_chain(accels_.acc(), chain_acc_, index());
    chain_advance(aux_vel_, chain_aux_vel_, chain_acc_, step_size);
  }

  template<typename Particles, typename Interactions>
  void ChainSystem<Particles, Interactions>::kick_real_vel(Scalar step_size) {
    std::swap(aux_vel_, ptcl_.vel());
    std::swap(chain_aux_vel_, chain_vel());
    interactions_.eval_extra_vel_dep_acc(*this, accels_.ext_vel_dep_acc());
    std::swap(aux_vel_, ptcl_.vel());
    std::swap(chain_aux_vel_, chain_vel());

    calc::coord_add(accels_.acc(), accels_.tot_vel_indep_acc(), accels_.ext_vel_dep_acc());
    Chain::calc_chain(accels_.acc(), chain_acc_, index());
    chain_advance(ptcl_.vel(), chain_vel_(), chain_acc_, step_size);
  }
}
#endif //SPACEHUB_ARCHAIN_HPP