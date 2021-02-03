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
#include "../interaction/interaction.hpp"
#include "../spacehub-concepts.hpp"

namespace space::particle_system {

    /*---------------------------------------------------------------------------*\
        Class SimpleSystem Declaration
    \*---------------------------------------------------------------------------*/
    /**
     * @tparam Particles
     * @tparam Interactions
     */
    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
    class SimpleSystem : public Particles {
       public:
        // Type members
        SPACEHUB_USING_TYPE_SYSTEM_OF(Particles);

        using Particle = typename Particles::Particle;

        using Interaction = Interactions;

        // Constructors
        SPACEHUB_MAKE_CONSTRUCTORS(SimpleSystem, delete, default, default, default, default);

        /**
         *
         * @tparam STL
         * @param time
         * @param particle_set
         */
        template <CONCEPT_PARTICLE_CONTAINER STL>
        SimpleSystem(Scalar time, STL const &particle_set);

        // Public methods

        /**
         * @param dt
         */
        void advance_time(Scalar dt);

        /**
         *
         * @param velocity
         * @param step_size
         */
        template <typename GenVectorArray>
        void advance_pos(Scalar step_size, GenVectorArray const &velocity);

        /**
         *
         * @param acceleration
         * @param step_size
         */
        template <typename GenVectorArray>
        void advance_vel(Scalar step_size, GenVectorArray const &acceleration);

        /**
         *
         * @param acceleration
         */
        template <typename GenVectorArray>
        void evaluate_acc(GenVectorArray &acceleration) const;

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

        void post_iter_process() {};

        /**
         *
         * @tparam STL
         * @param stl_ranges
         */
        template <typename ScalarIterable>
        void write_to_scalar_array(ScalarIterable &stl_ranges);

        /**
         *
         * @tparam STL
         * @param stl_ranges
         */
        template <typename ScalarIterable>
        void read_from_scalar_array(ScalarIterable const &stl_ranges);

        /**
         * @brief
         *
         * @tparam STL
         * @param stl_ranges
         */
        template <typename ScalarIterable>
        void evaluate_general_derivative(ScalarIterable &stl_ranges);

        size_t variable_number() const;

        // Friend functions
        template <CONCEPT_PARTICLES P, CONCEPT_INTERACTION F>
        friend std::ostream &operator<<(std::ostream &os, SimpleSystem<P, F> const &ps);

        template <CONCEPT_PARTICLES P, CONCEPT_INTERACTION F>
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
        // Particles ptcl_;

        interactions::InteractionData <Interactions, VectorArray> accels_;

        std::conditional_t<Interactions::ext_vel_dep, AdVectorArray, Empty> aux_vel_;
    };
}  // namespace space::particle_system

namespace space::particle_system {
    /*---------------------------------------------------------------------------*\
        Class SimpleSystem Implementation
    \*---------------------------------------------------------------------------*/
    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
    template <CONCEPT_PARTICLE_CONTAINER STL>
    SimpleSystem<Particles, Interactions>::SimpleSystem(Scalar
    time,
    const STL &particle_set
    )
    :
    Particles(time, particle_set
    ),
    accels_(particle_set
    .

    size()

    ) {
    if constexpr (Interactions::ext_vel_dep) {
    aux_vel_ = this->vel();
}
}

template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
template <typename ScalarIterable>
void SimpleSystem<Particles, Interactions>::read_from_scalar_array(const ScalarIterable &stl_ranges) {
    auto begin = stl_ranges.begin();
    this->time() = *begin;
    size_t len = this->number() * 3;
    auto pos_begin = begin + 1;
    auto pos_end = pos_begin + len;
    auto vel_begin = pos_end;
    auto vel_end = vel_begin + len;
    load_to_coords(pos_begin, pos_end, this->pos());
    load_to_coords(vel_begin, vel_end, this->vel());
    if constexpr (Interactions::ext_vel_dep) {
        auto aux_vel_begin = vel_end;
        auto aux_vel_end = aux_vel_begin + len;
        load_to_coords(aux_vel_begin, aux_vel_end, aux_vel_);
    }
}

template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
template <typename ScalarIterable>
void SimpleSystem<Particles, Interactions>::write_to_scalar_array(ScalarIterable &stl_ranges) {
    stl_ranges.clear();
    stl_ranges.reserve(this->number() * 3 * (2 + static_cast<size_t>(Interactions::ext_vel_dep)) + 1);
    stl_ranges.emplace_back(this->time());
    add_coords_to(stl_ranges, this->pos());
    add_coords_to(stl_ranges, this->vel());
    if constexpr (Interactions::ext_vel_dep) {
        add_coords_to(stl_ranges, aux_vel_);
    }
}

template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
template <typename ScalarIterable>
void SimpleSystem<Particles, Interactions>::evaluate_general_derivative(ScalarIterable &stl_ranges) {
    stl_ranges.clear();
    stl_ranges.reserve(this->number() * 3 * (2 + static_cast<size_t>(Interactions::ext_vel_dep)) + 1);
    stl_ranges.emplace_back(1);              // dt/dh
    add_coords_to(stl_ranges, this->vel());  // dp/dt
    Interactions::eval_acc(*this, this->accels_.acc());
    add_coords_to(stl_ranges, this->accels_.acc());  // dv/dt
    if constexpr (Interactions::ext_vel_dep) {
        add_coords_to(stl_ranges, this->accels_.acc());  // dw/dt
    }
}

template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
size_t SimpleSystem<Particles, Interactions>::variable_number() const {
    return this->number() * 3 * (2 + static_cast<size_t>(Interactions::ext_vel_dep)) + 1;
}

template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
void SimpleSystem<Particles, Interactions>::pre_iter_process() {
    if constexpr (Interactions::ext_vel_dep) {
        aux_vel_ = this->vel();
    }
}

template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
void SimpleSystem<Particles, Interactions>::kick(Scalar
step_size) {
if constexpr (Interactions::ext_vel_dep) {
Scalar half_step = 0.5 * step_size;

eval_vel_indep_acc();

kick_real_vel(half_step);
kick_pseu_vel(step_size);
kick_real_vel(half_step);
} else {
Interactions::eval_acc(*this, accels_.

acc()

);
calc::array_advance(this->

vel(), accels_

.

acc(), step_size

);
}
}

template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
void SimpleSystem<Particles, Interactions>::drift(Scalar
step_size) {
this->

time()

+=
step_size;
calc::array_advance(this->

pos(),

this->

vel(), step_size

);
}

template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
template <typename GenVectorArray>
void SimpleSystem<Particles, Interactions>::evaluate_acc(GenVectorArray &acceleration) const {
    Interactions::eval_acc(*this, acceleration);
}

template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
template <typename GenVectorArray>
void SimpleSystem<Particles, Interactions>::advance_vel(Scalar
step_size,
GenVectorArray const &acceleration
) {
calc::array_advance(this->

vel(), acceleration, step_size

);
}

template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
template <typename GenVectorArray>
void SimpleSystem<Particles, Interactions>::advance_pos(Scalar
step_size,
GenVectorArray const &velocity
) {
calc::array_advance(this->

pos(), velocity, step_size

);
}

template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
void SimpleSystem<Particles, Interactions>::advance_time(Scalar
dt) {
this->

time()

+=
dt;
}

template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
void SimpleSystem<Particles, Interactions>::kick_real_vel(Scalar
step_size) {
std::swap(aux_vel_,
this->

vel()

);
Interactions::eval_extra_vel_dep_acc(*this, accels_.

ext_vel_dep_acc()

);
std::swap(aux_vel_,
this->

vel()

);
calc::array_add(accels_
.

acc(), accels_

.

tot_vel_indep_acc(), accels_

.

ext_vel_dep_acc()

);
calc::array_advance(this->

vel(), accels_

.

acc(), step_size

);
}

template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
void SimpleSystem<Particles, Interactions>::kick_pseu_vel(Scalar
step_size) {
Interactions::eval_extra_vel_dep_acc(*this, accels_.

ext_vel_dep_acc()

);
calc::array_add(accels_
.

acc(), accels_

.

tot_vel_indep_acc(), accels_

.

ext_vel_dep_acc()

);
calc::array_advance(aux_vel_, accels_
.

acc(), step_size

);
}

template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
void SimpleSystem<Particles, Interactions>::eval_vel_indep_acc() {
    Interactions::eval_newtonian_acc(*this, accels_.tot_vel_indep_acc());
    if constexpr (Interactions::ext_vel_indep) {
        Interactions::eval_extra_vel_indep_acc(*this, accels_.ext_vel_indep_acc());
        calc::array_add(accels_.tot_vel_indep_acc(), accels_.tot_vel_indep_acc(), accels_.ext_vel_indep_acc());
    }
}

template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
std::ostream &operator<<(std::ostream &os, SimpleSystem <Particles, Interactions> const &ps) {
    os << static_cast<Particles>(ps);
    return os;
}

template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
std::istream &operator>>(std::istream &is, SimpleSystem <Particles, Interactions> &ps) {
    is >> static_cast<Particles>(ps);
    return is;
}

}  // namespace space::particle_system
