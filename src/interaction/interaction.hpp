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
    template <CONCEPT_FORCE InternalForce, CONCEPT_FORCE ExtraVelDepForce = void,
              CONCEPT_FORCE ExtraVelIndepForce = void>
    class Interactions {
       public:
        static constexpr bool ext_vel_dep{!std::is_same_v<ExtraVelDepForce, void>};

        static constexpr bool ext_vel_indep{!std::is_same_v<ExtraVelIndepForce, void>};

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

    template <typename Arg, typename... Args>
    struct ForceAdd {
        template <CONCEPT_PARTICLES_DATA Particles>
        static void add_acc_to(Particles const &particles, typename Particles::VectorArray &acceleration) {
            Arg::add_acc_to(particles, acceleration);
            if constexpr (sizeof...(Args) > 0) {
                ForceAdd<Args...>::add_acc_to(particles, acceleration);
            }
        }
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

    /*---------------------------------------------------------------------------*\
        Class Interactions Implementation
    \*---------------------------------------------------------------------------*/

    template <CONCEPT_FORCE InternalForce, CONCEPT_FORCE ExtraVelDepForce, CONCEPT_FORCE ExtraVelIndepForce>
    template <CONCEPT_PARTICLES_DATA Particles>
    void Interactions<InternalForce, ExtraVelDepForce, ExtraVelIndepForce>::eval_acc(
        const Particles &particles, typename Particles::VectorArray &acceleration) {
        calc::array_set_zero(acceleration);
        InternalForce::add_acc_to(particles, acceleration);
        if constexpr (ext_vel_dep) {
            ExtraVelDepForce::add_acc_to(particles, acceleration);
        }
        if constexpr (ext_vel_indep) {
            ExtraVelIndepForce::add_acc_to(particles, acceleration);
        }
    }

    template <CONCEPT_FORCE InternalForce, CONCEPT_FORCE ExtraVelDepForce, CONCEPT_FORCE ExtraVelIndepForce>
    template <CONCEPT_PARTICLES_DATA Particles>
    void Interactions<InternalForce, ExtraVelDepForce, ExtraVelIndepForce>::eval_extra_vel_dep_acc(
        const Particles &particles, typename Particles::VectorArray &acceleration) {
        if constexpr (ext_vel_dep) {
            calc::array_set_zero(acceleration);
            ExtraVelDepForce::add_acc_to(particles, acceleration);
        }
    }

    template <CONCEPT_FORCE InternalForce, CONCEPT_FORCE ExtraVelDepForce, CONCEPT_FORCE ExtraVelIndepForce>
    template <CONCEPT_PARTICLES_DATA Particles>
    void Interactions<InternalForce, ExtraVelDepForce, ExtraVelIndepForce>::eval_extra_vel_indep_acc(
        const Particles &particles, typename Particles::VectorArray &acceleration) {
        if constexpr (ext_vel_indep) {
            calc::array_set_zero(acceleration);
            ExtraVelIndepForce::add_acc_to(particles, acceleration);
        }
    }

    template <CONCEPT_FORCE InternalForce, CONCEPT_FORCE ExtraVelDepForce, CONCEPT_FORCE ExtraVelIndepForce>
    template <CONCEPT_PARTICLES_DATA Particles>
    void Interactions<InternalForce, ExtraVelDepForce, ExtraVelIndepForce>::eval_newtonian_acc(
        const Particles &particles, typename Particles::VectorArray &acceleration) {
        calc::array_set_zero(acceleration);
        InternalForce::add_acc_to(particles, acceleration);
    }

    template <typename Interactions, typename VectorArray>
    InteractionData<Interactions, VectorArray>::InteractionData(size_t size)
        : acc_{size}, newtonian_acc_{size}, tot_vel_indep_acc_{size} {
        if constexpr (Interactions::ext_vel_indep) {
            ext_vel_indep_acc_.resize(size);
        }
        if constexpr (Interactions::ext_vel_dep) {
            ext_vel_dep_acc_.resize(size);
        }
    }
}  // namespace space::interactions
