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

#define COMPILER_VERSION 201703L

#if __cplusplus > COMPILER_VERSION

#include <concepts>
#include <type_traits>

#include "type-class.hpp"
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
        ->std::convertible_to<typename T::Scalar>;
        { p.pos }
        ->std::convertible_to<typename T::Vector>;
        { p.vel }
        ->std::convertible_to<typename T::Vector>;
    };

    template <typename T>
    concept ParticleContainer = Container<T> &&Particle<typename T::value_type>;

    template <typename T>
    concept ParticlesData = requires(T p, size_t index) {
        typename T::Particle;
        typename T::Vector;
        typename T::Scalar;
        typename T::VectorArray;
        typename T::ScalarArray;
        typename T::IdxArray;

        requires Particle<typename T::Particle>;

        { p.time() }
        ->std::common_reference_with<typename T::Scalar &>;
        { p.idn() }
        ->std::common_reference_with<typename T::IdxArray &>;
        { p.mass() }
        ->std::common_reference_with<typename T::ScalarArray &>;
        { p.pos() }
        ->std::common_reference_with<typename T::VectorArray &>;
        { p.vel() }
        ->std::common_reference_with<typename T::VectorArray &>;
        { p.idn(index) }
        ->std::common_reference_with<size_t &>;
        { p.mass(index) }
        ->std::common_reference_with<typename T::Scalar &>;
        { p.pos(index) }
        ->std::common_reference_with<typename T::Vector &>;
        { p.vel(index) }
        ->std::common_reference_with<typename T::Vector &>;
        { p.number() }
        ->std::convertible_to<size_t> const;
    };

    template <typename T>
    concept ParticlesResizable = requires(T p, size_t index) {
        typename T::Particle;
        typename T::Vector;
        typename T::Scalar;
        typename T::VectorArray;
        typename T::ScalarArray;
        typename T::IdxArray;

        requires Particle<typename T::Particle>;
        { p.capacity() }
        ->std::convertible_to<size_t> const;
        { p.emplace_back(INSTANCE(typename T::Particle &)) }
        ->std::same_as<void>;
        { p.reserve(index) }
        ->std::same_as<void>;
        { p.resize(index) }
        ->std::same_as<void>;
        { p.clear() }
        ->std::same_as<void>;
    };

    template <typename T>
    concept Particles = ParticlesData<T> &&ParticlesResizable<T> &&requires(T p) {
        { p.column_names() }
        ->std::same_as<std::string>;
    };

    template <typename T>
    concept ParticleSystem = ParticlesData<T> &&requires(T p, size_t index) {
        typename T::Particle;
        typename T::Vector;
        typename T::Scalar;
        typename T::VectorArray;
        typename T::ScalarArray;
        typename T::IdxArray;
        requires Particles<std::remove_const_t<std::remove_reference_t<decltype(p.particles())>>>;

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
        { p.write_to_scalar_array(INSTANCE(typename T::ScalarArray &)) }
        ->std::same_as<void>;
        { p.read_from_scalar_array(INSTANCE(typename T::ScalarArray const &)) }
        ->std::same_as<void>;
        { p.column_names() }
        ->std::same_as<std::string>;
    };

    struct TestParticle {
        SPACEHUB_USING_TYPE_SYSTEM_OF(Types<double>);
        Scalar mass;
        Vector pos;
        Vector vel;
    };
    struct TestParticles {
        SPACEHUB_USING_TYPE_SYSTEM_OF(TestParticle)
        using Particle = TestParticle;
        Scalar &time();
        IdxArray &idn();
        ScalarArray &mass();
        VectorArray &pos();
        VectorArray &vel();
        size_t &idn(size_t);
        Scalar &mass(size_t);
        Vector &pos(size_t);
        Vector &vel(size_t);
        size_t number() const;
        size_t capacity() const;
        void emplace_back(Particle &);
        void reserve(size_t);
        void resize(size_t);
        void clear();
    };

    template <typename T>
    concept NewtonianInteraction = requires(TestParticles p) {
        { T::eval_acc(p, INSTANCE(typename TestParticles::VectorArray &)) }
        ->std::same_as<void>;

        { T::eval_newtonian_acc(p, INSTANCE(typename TestParticles::VectorArray &)) }
        ->std::same_as<void>;
    };

    template <typename T>
    concept ExtraVelDepInteraction = requires(TestParticles p) {
        // requires Particles<typename T::Particles>;

        { T::eval_extra_vel_dep_acc(p, INSTANCE(typename TestParticles::VectorArray &)) }
        ->std::same_as<void>;
    };

    template <typename T>
    concept ExtraVelIndepInteraction = requires(TestParticles p) {
        // requires Particles<typename T::Particles>;

        { T::eval_extra_vel_indep_acc(p, INSTANCE(typename TestParticles::VectorArray &)) }
        ->std::same_as<void>;
    };

    template <typename T>
    concept Interaction = NewtonianInteraction<T> &&
                          (!T::ext_vel_dep || (T::ext_vel_dep && ExtraVelDepInteraction<T>)) &&
                          (!T::ext_vel_indep || (T::ext_vel_indep && ExtraVelIndepInteraction<T>));

    template <typename T>
    concept Force = std::same_as<T, void> || requires(TestParticles p) {
        { T::vel_dependent }
        ->std::convertible_to<bool>;
        { T::add_acc_to(p, INSTANCE(typename TestParticles::VectorArray &)) }
        ->std::same_as<void>;
    };

    template <typename T>
    concept StepControler = requires(T c, std::tuple<typename T::Scalar, typename T::Scalar> tup) {
        typename T::Scalar;
        { c.next_step_size(INSTANCE(typename T::Scalar), INSTANCE(typename T::Scalar), INSTANCE(typename T::Scalar)) }
        ->std::convertible_to<typename T::Scalar>;
        { c.next_step_size(INSTANCE(typename T::Scalar), INSTANCE(typename T::Scalar), tup) }
        ->std::convertible_to<typename T::Scalar>;
    };
    // template <typename T>
    // concept Integrator = requires(T i){};

}  // namespace space::concepts

#define CONCEPT_PARTICLE concepts::Particle
#define CONCEPT_PARTICLES concepts::Particles
#define CONCEPT_PARTICLE_SYSTEM concepts::ParticleSystem
#define CONCEPT_CONTAINER concepts::Container
#define CONCEPT_PARTICLES_DATA concepts::ParticlesData
#define CONCEPT_INTERACTION concepts::Interaction
#define CONCEPT_PARTICLE_CONTAINER concepts::ParticleContainer
#define CONCEPT_FORCE concepts::Force
#define CONCEPT_STEP_CONTROLER concepts::StepControler

#else

#define CONCEPT_PARTICLE typename
#define CONCEPT_PARTICLES typename
#define CONCEPT_PARTICLE_SYSTEM typename
#define CONCEPT_CONTAINER typename
#define CONCEPT_PARTICLES_DATA typename
#define CONCEPT_INTERACTION typename
#define CONCEPT_PARTICLE_CONTAINER typename
#define CONCEPT_FORCE typename
#define CONCEPT_STEP_CONTROLER typename
#endif
