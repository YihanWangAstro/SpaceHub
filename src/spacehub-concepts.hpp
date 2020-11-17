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
 * @file spacehub-concepts.hpp
 *
 * Header file.
 */
#pragma once

#include <concepts>
namespace space::concepts {

#define INSTANCE(X) std::declval<X>()

template <typename T>
concept Container = requires(T t) {
  typename T::value_type;
  { t.size() }
  ->std::convertible_to<size_t>;
  {t.begin()};
  {t.end()};
};

template <typename T>
concept Particle = requires(T p) {
  typename T::Scalar;
  typename T::Vector;

  { p.mass }
  ->std::convertible_to<typename T::Scalar &>;
  { p.pos }
  ->std::convertible_to<typename T::Vector &>;
  { p.vel }
  ->std::convertible_to<typename T::Vector &>;
};

template <typename T>
concept ParticleContainer = Container<T> &&Particle<typename T::value_type>;

template <typename T>
concept Particles = requires(T p, size_t index) {
  typename T::Particle;
  typename T::Vector;
  typename T::Scalar;
  typename T::VectorArray;
  typename T::ScalarArray;
  typename T::IdxArray;

  { p.time() }
  ->std::convertible_to<typename T::Scalar &>;
  { p.idn() }
  ->std::convertible_to<typename T::IdxArray &>;
  { p.mass() }
  ->std::convertible_to<typename T::ScalarArray &>;
  { p.pos() }
  ->std::convertible_to<typename T::VectorArray &>;
  { p.vel() }
  ->std::convertible_to<typename T::VectorArray &>;
  { p.idn(index) }
  ->std::convertible_to<size_t &>;
  { p.mass(index) }
  ->std::convertible_to<typename T::Scalar &>;
  { p.pos(index) }
  ->std::convertible_to<typename T::Vector &>;
  { p.vel(index) }
  ->std::convertible_to<typename T::Vector &>;
  { p.number() }
  ->std::convertible_to<size_t> const;
  { p.capacity() }
  ->std::convertible_to<size_t> const;
  //{ p.emplace_back(declval<typename T::Particle>()) }
  { p.emplace_back(INSTANCE(typename T::Particle)) }
  ->std::same_as<void>;
  { p.reserve(index) }
  ->std::same_as<void>;
  { p.resize(index) }
  ->std::same_as<void>;
  { p.clear() }
  ->std::same_as<void>;
};

template <typename T>
concept ParticleSystem = requires(T p, size_t index) {
  typename T::Particle;
  typename T::Vector;
  typename T::Scalar;
  typename T::VectorArray;
  typename T::ScalarArray;
  typename T::IdxArray;

  requires Particles<decltype(p.particles())>;

  { p.advance_time(INSTANCE(typename T::Scalar)) }
  ->std::same_as<void>;
  { p.advance_pos(INSTANCE(typename T::Scalar), INSTANCE(typename T::VectorArray const &)) }
  ->std::same_as<void>;
  { p.advance_vel(INSTANCE(typename T::Scalar), INSTANCE(typename T::VectorArray const &)) }
  ->std::same_as<void>;
  { p.evaluate_acc(INSTANCE(typename T::VectorArray &)) }
  ->std::same_as<void>;
  { p.drift(INSTANCE(typename T::Scalar)) }
  ->std::same_as<void>;
  { p.kick(INSTANCE(typename T::Scalar)) }
  ->std::same_as<void>;
  { p.pre_iter_process() }
  ->std::same_as<void>;
  { p.post_iter_process() }
  ->std::same_as<void>;
  { p.to_linear_container(INSTANCE(typename T::ScalarArray &)) }
  ->std::same_as<void>;
  { p.load_from_linear_container(INSTANCE(typename T::ScalarArray const &)) }
  ->std::same_as<void>;
};

template <typename T>
concept Interactions = requires(T i) {
  requires Particles<typename T::Particles>;

  { i.eval_acc(INSTANCE(typename T::Particles const &), INSTANCE(typename T::Particles::VectorArray &)) }
  ->std::same_as<void>;

  { i.eval_extra_vel_dep_acc(INSTANCE(typename T::Particles const &), INSTANCE(typename T::Particles::VectorArray &)) }
  ->std::same_as<void>;

  {
    i.eval_extra_vel_indep_acc(INSTANCE(typename T::Particles const &), INSTANCE(typename T::Particles::VectorArray &))
  }
  ->std::same_as<void>;

  { i.eval_newtonian_acc(INSTANCE(typename T::Particles const &), INSTANCE(typename T::Particles::VectorArray &)) }
  ->std::same_as<void>;
};

}  // namespace space::concepts
