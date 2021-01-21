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
#pragma once

#include <type_traits>

#include "../core-computation.hpp"
#include "chain.hpp"

namespace space::particle_system {

    /*---------------------------------------------------------------------------*\
        Class ChainSystem Declaration
    \*---------------------------------------------------------------------------*/
    /**
     *
     * @tparam Particles
     * @tparam Interactions
     */
    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
    class ChainSystem : public Particles {
       public:
        // Type members
        SPACEHUB_USING_TYPE_SYSTEM_OF(Particles);

        using Particle = typename Particles::Particle;

        // Constructors
        SPACEHUB_MAKE_CONSTRUCTORS(ChainSystem, delete, default, default, default, default);

        template <CONCEPT_PARTICLE_CONTAINER STL>
        ChainSystem(Scalar time, STL const &particle_set);

        // Public methods
        SPACEHUB_STD_ACCESSOR(VectorArray, chain_pos, chain_pos_);

        SPACEHUB_STD_ACCESSOR(VectorArray, chain_vel, chain_vel_);

        SPACEHUB_STD_ACCESSOR(IdxArray, index, index_);

        void advance_time(Scalar dt);

        void advance_pos(Scalar step_size, VectorArray const &velocity);

        void advance_vel(Scalar step_size, VectorArray const &acceleration);

        void evaluate_acc(VectorArray &acceleration) const;

        void drift(Scalar step_size);

        void kick(Scalar step_size);

        void pre_iter_process();

        void post_iter_process();

        template <typename STL>
        void write_to_scalar_array(STL &stl_ranges);

        template <typename STL>
        void read_from_scalar_array(STL const &stl_ranges);

       private:
        // Private methods
        void chain_advance(VectorArray &var, VectorArray &chain_var, VectorArray &chain_increment, Scalar step_size);

        void eval_vel_indep_acc();

        void kick_pseu_vel(Scalar step_size);

        void kick_real_vel(Scalar step_size);

        // Friend functions
        template <CONCEPT_PARTICLES P, CONCEPT_INTERACTION F>
        friend std::ostream &operator<<(std::ostream &os, ChainSystem<P, F> const &ps);

        template <CONCEPT_PARTICLES P, CONCEPT_INTERACTION F>
        friend std::istream &operator>>(std::istream &is, ChainSystem<P, F> &ps);

       private:
        // Private members
        interactions::InteractionData<Interactions, VectorArray> accels_;
        VectorArray chain_pos_;
        VectorArray chain_vel_;
        VectorArray chain_acc_;

        IdxArray index_;
        IdxArray new_index_;

        std::conditional_t<Interactions::ext_vel_dep, VectorArray, Empty> aux_vel_;
        std::conditional_t<Interactions::ext_vel_dep, VectorArray, Empty> chain_aux_vel_;
    };

    /*---------------------------------------------------------------------------*\
        Class ChainSystem Implementation
    \*---------------------------------------------------------------------------*/
    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
    template <CONCEPT_PARTICLE_CONTAINER STL>
    ChainSystem<Particles, Interactions>::ChainSystem(Scalar time, const STL &particle_set)
        : Particles(time, particle_set),
          chain_pos_(particle_set.size()),
          chain_vel_(particle_set.size()),
          index_(particle_set.size()),
          new_index_(particle_set.size()),
          accels_(particle_set.size()),
          chain_acc_(particle_set.size()) {
        Chain::calc_chain_index(this->pos(), index_);
        Chain::calc_chain(this->pos(), chain_pos(), index_);
        Chain::calc_chain(this->vel(), chain_vel(), index_);

        if constexpr (Interactions::ext_vel_dep) {
            aux_vel_ = this->vel();
            chain_aux_vel_ = chain_vel_;
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
    void ChainSystem<Particles, Interactions>::advance_time(Scalar dt) {
        this->time() += dt;
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
    void ChainSystem<Particles, Interactions>::advance_pos(Scalar step_size, VectorArray const &velocity) {
        Chain::calc_chain(velocity, chain_vel(), index());
        chain_advance(this->pos(), chain_pos(), chain_vel(), step_size);
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
    void ChainSystem<Particles, Interactions>::advance_vel(Scalar step_size, VectorArray const &acceleration) {
        Chain::calc_chain(acceleration, chain_acc_, index());
        chain_advance(this->vel(), chain_vel(), chain_acc_, step_size);
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
    void ChainSystem<Particles, Interactions>::evaluate_acc(VectorArray &acceleration) const {
        Interactions::eval_acc(*this, acceleration);
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
    void ChainSystem<Particles, Interactions>::drift(Scalar step_size) {
        this->time() += step_size;
        chain_advance(this->pos(), chain_pos(), chain_vel(), step_size);
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
    void ChainSystem<Particles, Interactions>::kick(Scalar step_size) {
        if constexpr (Interactions::ext_vel_dep) {
            Scalar half_step = 0.5 * step_size;
            eval_vel_indep_acc();
            kick_real_vel(half_step);
            kick_pseu_vel(step_size);
            kick_real_vel(half_step);
        } else {
            Interactions::eval_acc(*this, accels_.acc());
            advance_vel(step_size, accels_.acc());
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
    void ChainSystem<Particles, Interactions>::pre_iter_process() {
        if constexpr (Interactions::ext_vel_dep) {
            aux_vel_ = this->vel();
            chain_aux_vel_ = chain_vel_;
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
    void ChainSystem<Particles, Interactions>::post_iter_process() {
        Chain::calc_chain_index(this->pos(), new_index_);
        if (new_index_ != index_) {
            Chain::update_chain(chain_pos_, this->pos(), index_, new_index_);
            Chain::calc_cartesian(this->mass(), chain_pos_, this->pos(), new_index_);
            Chain::update_chain(chain_vel_, this->vel(), index_, new_index_);
            Chain::calc_cartesian(this->mass(), chain_vel_, this->vel(), new_index_);
            index_ = new_index_;
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
    template <typename STL>
    void ChainSystem<Particles, Interactions>::write_to_scalar_array(STL &stl_ranges) {
        stl_ranges.clear();
        stl_ranges.reserve(this->number() * 3 * (2 + static_cast<size_t>(Interactions::ext_vel_dep)) + 1);
        stl_ranges.emplace_back(this->time());
        add_coords_to(stl_ranges, chain_pos_);
        add_coords_to(stl_ranges, chain_vel_);
        if constexpr (Interactions::ext_vel_dep) {
            add_coords_to(stl_ranges, chain_aux_vel_);
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
    template <typename STL>
    void ChainSystem<Particles, Interactions>::read_from_scalar_array(const STL &stl_ranges) {
        auto begin = stl_ranges.begin();
        this->time() = *begin;
        size_t len = this->number() * 3;
        auto pos_begin = begin + 1;
        auto pos_end = pos_begin + len;
        auto vel_begin = pos_end;
        auto vel_end = vel_begin + len;

        load_to_coords(pos_begin, pos_end, chain_pos_);
        load_to_coords(vel_begin, vel_end, chain_vel_);

        Chain::calc_cartesian(this->mass(), chain_pos_, this->pos(), index_);
        Chain::calc_cartesian(this->mass(), chain_vel_, this->vel(), index_);
        if constexpr (Interactions::ext_vel_dep) {
            auto aux_vel_begin = vel_end;
            auto aux_vel_end = aux_vel_begin + len;
            load_to_coords(aux_vel_begin, aux_vel_end, chain_aux_vel_);
            Chain::calc_cartesian(this->mass(), chain_aux_vel_, aux_vel_, index_);
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
    std::istream &operator>>(std::istream &is, ChainSystem<Particles, Interactions> &ps) {
        is >> static_cast<Particles>(ps);
        return is;
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
    std::ostream &operator<<(std::ostream &os, const ChainSystem<Particles, Interactions> &ps) {
        os << static_cast<Particles>(ps);
        return os;
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
    void ChainSystem<Particles, Interactions>::chain_advance(VectorArray &var, VectorArray &chain_var,
                                                             VectorArray &chain_increment, Scalar step_size) {
        calc::array_advance(chain_var, chain_increment, step_size);
        Chain::calc_cartesian(this->mass(), chain_var, var, index());
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
    void ChainSystem<Particles, Interactions>::eval_vel_indep_acc() {
        Interactions::eval_newtonian_acc(*this, accels_.tot_vel_indep_acc());
        if constexpr (Interactions::ext_vel_indep) {
            Interactions::eval_extra_vel_indep_acc(*this, accels_.ext_vel_indep_acc());
            calc::array_add(accels_.tot_vel_indep_acc(), accels_.tot_vel_indep_acc(), accels_.ext_vel_indep_acc());
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
    void ChainSystem<Particles, Interactions>::kick_pseu_vel(Scalar step_size) {
        Interactions::eval_extra_vel_dep_acc(*this, accels_.ext_vel_dep_acc());
        calc::array_add(accels_.acc(), accels_.tot_vel_indep_acc(), accels_.ext_vel_dep_acc());
        Chain::calc_chain(accels_.acc(), chain_acc_, index());
        chain_advance(aux_vel_, chain_aux_vel_, chain_acc_, step_size);
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
    void ChainSystem<Particles, Interactions>::kick_real_vel(Scalar step_size) {
        std::swap(aux_vel_, this->vel());
        std::swap(chain_aux_vel_, chain_vel());
        Interactions::eval_extra_vel_dep_acc(*this, accels_.ext_vel_dep_acc());
        std::swap(aux_vel_, this->vel());
        std::swap(chain_aux_vel_, chain_vel());

        calc::array_add(accels_.acc(), accels_.tot_vel_indep_acc(), accels_.ext_vel_dep_acc());
        Chain::calc_chain(accels_.acc(), chain_acc_, index());
        chain_advance(this->vel(), chain_vel(), chain_acc_, step_size);
    }
}  // namespace space::particle_system
