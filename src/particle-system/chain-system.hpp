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
    the terms of the GPL-3.0 License. SpaceHub is distributed in the hope that it
    will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GPL-3.0 License
    for more details. You should have received a copy of the GPL-3.0 License along
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
#include "../type-class.hpp"
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

        using Interaction = Interactions;

        static constexpr bool ext_vel_dep{Interactions::ext_vel_dep};

        static constexpr bool ext_vel_indep{Interactions::ext_vel_indep};

        // Constructors
        SPACEHUB_MAKE_CONSTRUCTORS(ChainSystem, delete, default, default, default, default);

        template <CONCEPT_PARTICLE_CONTAINER STL>
        ChainSystem(Scalar time, STL const &particle_set);

        // Public methods
        SPACEHUB_ARRAY_ACCESSOR(StateVectorArray, chain_pos, chain_pos_);

        SPACEHUB_ARRAY_ACCESSOR(StateVectorArray, chain_vel, chain_vel_);

        SPACEHUB_ARRAY_ACCESSOR(IdxArray, index, index_);

        SPACEHUB_ARRAY_ACCESSOR(StateScalarArray, increment, increment_);

        Scalar step_scale() const { return 1.0; };

        template <typename GenVectorArray>
        void evaluate_acc(GenVectorArray &acceleration) const;

        void drift(Scalar step_size);

        void kick(Scalar step_size);

        void pre_iter_process();

        void post_iter_process();

        template <typename ScalarIterable>
        void write_to_scalar_array(ScalarIterable &stl_ranges);

        template <typename ScalarIterable>
        void read_from_scalar_array(ScalarIterable const &stl_ranges);

        template <typename ScalarIterable>
        void evaluate_general_derivative(ScalarIterable &stl_ranges);

        inline void collect_increment(bool sync) { sync_increment_ = sync; };

        void clear_increment() { calc::array_set_zero(increment_); };

        [[nodiscard]] size_t variable_number() const;

        inline constexpr size_t time_offset() const { return 0; };
        inline constexpr size_t pos_offset() const { return 1; };
        inline constexpr size_t vel_offset() const { return this->number() * 3 + 1; };
        inline constexpr size_t auxi_vel_offset() const { return this->number() * 6 + 1; };

       private:
        // Private methods
        template <typename Array1, typename Array2, typename Array3>
        void chain_advance(Array1 &var, Array2 &chain_var, Array3 &chain_increment, Scalar step_size);

        void eval_vel_indep_acc();

        void kick_pseu_vel(Scalar step_size);

        void kick_real_vel(Scalar step_size);

        template <typename Array>
        void sync_pos_increment(Array const &inc, Scalar step_size);

        template <typename Array>
        void sync_vel_increment(Array const &inc, Scalar step_size);

        template <typename Array>
        void sync_auxi_vel_increment(Array const &inc, Scalar step_size);

        void sync_time_increment(Scalar phy_time);

        // Friend functions
        template <CONCEPT_PARTICLES P, CONCEPT_INTERACTION F>
        friend std::ostream &operator<<(std::ostream &os, ChainSystem<P, F> const &ps);

        template <CONCEPT_PARTICLES P, CONCEPT_INTERACTION F>
        friend std::istream &operator>>(std::istream &is, ChainSystem<P, F> &ps);

       private:
        // Private members
        force::InteractionData<Interactions, VectorArray> accels_;
        StateVectorArray chain_pos_;
        StateVectorArray chain_vel_;
        VectorArray chain_acc_;
        StateScalarArray increment_;
        IdxArray index_;
        IdxArray new_index_;

        std::conditional_t<Interactions::ext_vel_dep, StateVectorArray, Empty> aux_vel_;
        std::conditional_t<Interactions::ext_vel_dep, StateVectorArray, Empty> chain_aux_vel_;

        bool sync_increment_{false};
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
          chain_acc_(particle_set.size()),
          increment_(this->variable_number()) {
        Chain::calc_chain_index(this->pos(), index_);
        Chain::calc_chain(this->pos(), chain_pos(), index_);
        Chain::calc_chain(this->vel(), chain_vel(), index_);

        if constexpr (Interactions::ext_vel_dep) {
            aux_vel_ = this->vel();
            chain_aux_vel_ = chain_vel_;
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
    template <typename GenVectorArray>
    void ChainSystem<Particles, Interactions>::evaluate_acc(GenVectorArray &acceleration) const {
        Interactions::eval_acc(*this, acceleration);
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
    void ChainSystem<Particles, Interactions>::drift(Scalar step_size) {
        this->time() += step_size;
        chain_advance(this->pos(), chain_pos(), chain_vel(), step_size);
        sync_time_increment(step_size);
        sync_pos_increment(chain_vel(), step_size);
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
            Chain::calc_chain(accels_.acc(), chain_acc_, index());
            chain_advance(this->vel(), chain_vel(), chain_acc_, step_size);
            sync_vel_increment(chain_acc_, step_size);
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
    template <typename ScalarIterable>
    void ChainSystem<Particles, Interactions>::write_to_scalar_array(ScalarIterable &stl_ranges) {
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
    template <typename ScalarIterable>
    void ChainSystem<Particles, Interactions>::evaluate_general_derivative(ScalarIterable &stl_ranges) {
        stl_ranges.clear();
        stl_ranges.reserve(this->number() * 3 * (2 + static_cast<size_t>(Interactions::ext_vel_dep)) + 1);
        stl_ranges.emplace_back(1);             // dt/dt
        add_coords_to(stl_ranges, chain_vel_);  // dX/dt
        evaluate_acc(this->accels_.acc());
        Chain::calc_chain(this->accels_.acc(), chain_acc_, index());
        add_coords_to(stl_ranges, chain_acc_);  // dV/dt
        if constexpr (Interactions::ext_vel_dep) {
            add_coords_to(stl_ranges, chain_acc_);
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
    template <typename ScalarIterable>
    void ChainSystem<Particles, Interactions>::read_from_scalar_array(const ScalarIterable &stl_ranges) {
        if (stl_ranges.size() == this->variable_number()) {
            auto begin = stl_ranges.begin();
            this->time() = *(begin + time_offset());
            auto pos_begin = begin + pos_offset();
            auto pos_end = begin + vel_offset();
            auto vel_begin = begin + vel_offset();
            auto vel_end = begin + auxi_vel_offset();

            load_to_coords(pos_begin, pos_end, chain_pos_);
            load_to_coords(vel_begin, vel_end, chain_vel_);

            Chain::calc_cartesian(this->mass(), chain_pos_, this->pos(), index_);
            Chain::calc_cartesian(this->mass(), chain_vel_, this->vel(), index_);
            if constexpr (Interactions::ext_vel_dep) {
                auto aux_vel_begin = begin + auxi_vel_offset();
                auto aux_vel_end = stl_ranges.end();
                load_to_coords(aux_vel_begin, aux_vel_end, chain_aux_vel_);
                Chain::calc_cartesian(this->mass(), chain_aux_vel_, aux_vel_, index_);
            }
        } else {
            spacehub_abort("Wrong input array size!");
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>

    size_t ChainSystem<Particles, Interactions>::variable_number() const {
        return this->number() * 3 * (2 + static_cast<size_t>(Interactions::ext_vel_dep)) + 1;
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
    template <typename Array1, typename Array2, typename Array3>
    void ChainSystem<Particles, Interactions>::chain_advance(Array1 &var, Array2 &chain_var, Array3 &chain_increment,
                                                             Scalar step_size) {
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
        sync_auxi_vel_increment(chain_acc_, step_size);
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
        sync_vel_increment(chain_acc_, step_size);
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
    template <typename Array>
    void ChainSystem<Particles, Interactions>::sync_pos_increment(Array const &inc, Scalar step_size) {
        if (sync_increment_) {
            advance_scaled_coords_to(increment_.begin() + pos_offset(), inc, step_size);
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
    template <typename Array>
    void ChainSystem<Particles, Interactions>::sync_vel_increment(Array const &inc, Scalar step_size) {
        if (sync_increment_) {
            advance_scaled_coords_to(increment_.begin() + vel_offset(), inc, step_size);
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
    template <typename Array>
    void ChainSystem<Particles, Interactions>::sync_auxi_vel_increment(Array const &inc, Scalar step_size) {
        if (sync_increment_) {
            advance_scaled_coords_to(increment_.begin() + auxi_vel_offset(), inc, step_size);
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
    void ChainSystem<Particles, Interactions>::sync_time_increment(Scalar phy_time) {
        if (sync_increment_) {
            increment_[time_offset()] += phy_time;
        }
    }
}  // namespace space::particle_system
