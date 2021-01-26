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
 * @file interaction.hpp
 *
 * Header file.
 */
#pragma once

#include "../core-computation.hpp"
#include "../spacehub-concepts.hpp"
namespace space::interactions {

    /*---------------------------------------------------------------------------*\
        Class Interactions Declaration
    \*---------------------------------------------------------------------------*/
    template <CONCEPT_FORCE InternalForce, CONCEPT_FORCE... ExtraForce>
    class Interactions {
       public:
        static constexpr bool ext_vel_dep{(... || ExtraForce::vel_dependent)};

        static constexpr bool ext_vel_indep{(... || !ExtraForce::vel_dependent)};

        /**
         * Evaluate the total acceleration of the current state of a given particle system.
         *
         * @tparam Particles Type of the particle system.
         *
         * @param[in] particles The particle system need to be evaluated.
         * @param[out] acceleration The output of the evaluated acceleration.
         */
        template <CONCEPT_PARTICLES_DATA Particles>
        static void eval_acc(Particles const &particles, typename Particles::VectorArray &acceleration);

        /**
         * Evaluate the external acceleration of the current state of a given particle system.
         *
         * @brief
         *
         * @tparam Particles
         * @param particles
         * @param acceleration
         */
        template <CONCEPT_PARTICLES_DATA Particles>
        static void eval_extra_acc(Particles const &particles, typename Particles::VectorArray &acceleration);

        /**
         * Evaluate the external velocity dependent acceleration of the current state of a given particle system.
         *
         * @tparam Particles Type of the particle system.
         *
         * @param[in] particles The particle system need to be evaluated.
         * @param[out] acceleration The output of the evaluated acceleration.
         */
        template <CONCEPT_PARTICLES_DATA Particles>
        static void eval_extra_vel_dep_acc(Particles const &particles, typename Particles::VectorArray &acceleration);

        /**
         * Evaluate the external velocity independent acceleration of the current state of a given particle system.
         *
         * @tparam Particles Type of the particle system.
         *
         * @param[in] particles The particle system need to be evaluated.
         * @param[out] acceleration The output of the evaluated acceleration.
         */
        template <CONCEPT_PARTICLES_DATA Particles>
        static void eval_extra_vel_indep_acc(Particles const &particles, typename Particles::VectorArray &acceleration);

        /**
         * Evaluate the internal newtonian acceleration of the current state of a given particle system.
         *
         * @tparam Particles Type of the particle system.
         *
         * @param[in] particles The particle system need to be evaluated.
         * @param[out] acceleration The output of the evaluated acceleration.
         */
        template <CONCEPT_PARTICLES_DATA Particles>
        static void eval_newtonian_acc(Particles const &particles, typename Particles::VectorArray &acceleration);
    };

    template <typename Interactions, typename VectorArray>
    class InteractionData {
       public:
        // Constructors
        SPACEHUB_MAKE_CONSTRUCTORS(InteractionData, default, default, default, default, default);

        explicit InteractionData(size_t size);

        // Public methods

        SPACEHUB_ARRAY_ACCESSOR(VectorArray, acc, acc_);

        SPACEHUB_ARRAY_ACCESSOR(VectorArray, newtonian_acc, newtonian_acc_);

        SPACEHUB_ARRAY_ACCESSOR(VectorArray, tot_vel_indep_acc, tot_vel_indep_acc_);

        SPACEHUB_ARRAY_ACCESSOR(VectorArray, ext_vel_indep_acc, ext_vel_indep_acc_);

        SPACEHUB_ARRAY_ACCESSOR(VectorArray, ext_vel_dep_acc, ext_vel_dep_acc_);

       private:
        VectorArray acc_;

        VectorArray newtonian_acc_;

        VectorArray tot_vel_indep_acc_;

        std::conditional_t<Interactions::ext_vel_indep, VectorArray, Empty> ext_vel_indep_acc_;

        std::conditional_t<Interactions::ext_vel_dep, VectorArray, Empty> ext_vel_dep_acc_;
    };
    template <typename Arg, typename... Args>
    struct InvokeVelDepForce {
        template <CONCEPT_PARTICLES_DATA Particles>
        static void add_acc_to(Particles const &particles, typename Particles::VectorArray &acceleration) {
            if constexpr (Arg::vel_dependent) {
                Arg::add_acc_to(particles, acceleration);
            }
            if constexpr (sizeof...(Args) > 0) {
                InvokeVelDepForce<Args...>::add_acc_to(particles, acceleration);
            }
        }
    };

    template <typename Arg, typename... Args>
    struct InvokeVelIndepForce {
        template <CONCEPT_PARTICLES_DATA Particles>
        static void add_acc_to(Particles const &particles, typename Particles::VectorArray &acceleration) {
            if constexpr (!Arg::vel_dependent) {
                Arg::add_acc_to(particles, acceleration);
            }
            if constexpr (sizeof...(Args) > 0) {
                InvokeVelIndepForce<Args...>::add_acc_to(particles, acceleration);
            }
        }
    };
    /*---------------------------------------------------------------------------*\
        Class Interactions Implementation
    \*---------------------------------------------------------------------------*/

    template <CONCEPT_FORCE InternalForce, CONCEPT_FORCE... ExtraForce>
    template <CONCEPT_PARTICLES_DATA Particles>
    void Interactions<InternalForce, ExtraForce...>::eval_acc(const Particles &particles,
                                                              typename Particles::VectorArray &acceleration) {
        calc::array_set_zero(acceleration);
        InternalForce::add_acc_to(particles, acceleration);
        (ExtraForce::add_acc_to(particles, acceleration), ...);
    }

    template <CONCEPT_FORCE InternalForce, CONCEPT_FORCE... ExtraForce>
    template <CONCEPT_PARTICLES_DATA Particles>
    void Interactions<InternalForce, ExtraForce...>::eval_extra_acc(const Particles &particles,
                                                                    typename Particles::VectorArray &acceleration) {
        if constexpr (ext_vel_dep) {
            calc::array_set_zero(acceleration);
            InvokeVelDepForce<ExtraForce...>::add_acc_to(particles, acceleration);
        }
        if constexpr (ext_vel_indep) {
            InvokeVelIndepForce<ExtraForce...>::add_acc_to(particles, acceleration);
        }
    }

    template <CONCEPT_FORCE InternalForce, CONCEPT_FORCE... ExtraForce>
    template <CONCEPT_PARTICLES_DATA Particles>
    void Interactions<InternalForce, ExtraForce...>::eval_extra_vel_dep_acc(
        const Particles &particles, typename Particles::VectorArray &acceleration) {
        if constexpr (ext_vel_dep) {
            calc::array_set_zero(acceleration);
            InvokeVelDepForce<ExtraForce...>::add_acc_to(particles, acceleration);
        }
    }

    template <CONCEPT_FORCE InternalForce, CONCEPT_FORCE... ExtraForce>
    template <CONCEPT_PARTICLES_DATA Particles>
    void Interactions<InternalForce, ExtraForce...>::eval_extra_vel_indep_acc(
        const Particles &particles, typename Particles::VectorArray &acceleration) {
        if constexpr (ext_vel_indep) {
            calc::array_set_zero(acceleration);
            InvokeVelIndepForce<ExtraForce...>::add_acc_to(particles, acceleration);
        }
    }

    template <CONCEPT_FORCE InternalForce, CONCEPT_FORCE... ExtraForce>
    template <CONCEPT_PARTICLES_DATA Particles>
    void Interactions<InternalForce, ExtraForce...>::eval_newtonian_acc(const Particles &particles,
                                                                        typename Particles::VectorArray &acceleration) {
        calc::array_set_zero(acceleration);
        InternalForce::add_acc_to(particles, acceleration);
    }

    /*---------------------------------------------------------------------------*\
            Class InteractionData Implementation
    \*---------------------------------------------------------------------------*/
    template <typename Interactions, typename VectorArray>
    InteractionData<Interactions, VectorArray>::InteractionData(size_t size)
        : acc_(size), newtonian_acc_(size), tot_vel_indep_acc_(size) {
        if constexpr (Interactions::ext_vel_indep) {
            ext_vel_indep_acc_.resize(size);
        }
        if constexpr (Interactions::ext_vel_dep) {
            ext_vel_dep_acc_.resize(size);
        }
    }
}  // namespace space::interactions
