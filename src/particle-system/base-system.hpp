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
 * @file base-system.hpp
 *
 * Header file.
 */
#pragma once

#include <type_traits>

#include "../core-computation.hpp"
#include "../interaction/interaction.hpp"
#include "../spacehub-concepts.hpp"
#include "../type-class.hpp"
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

        static constexpr bool ext_vel_dep{Interactions::ext_vel_dep};

        static constexpr bool ext_vel_indep{Interactions::ext_vel_indep};

        // Constructors
        SPACEHUB_MAKE_CONSTRUCTORS(SimpleSystem, delete, default, default, default, default);

        SPACEHUB_ARRAY_ACCESSOR(StateScalarArray, increment, increment_);

        Scalar step_scale() const { return 1.0; };

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

        void post_iter_process(){};

        /**
         * @brief
         *
         * @tparam ScalarIterable
         * @param stl_ranges
         */
        template <typename ScalarIterable>
        void write_to_scalar_array(ScalarIterable &stl_ranges);

        /**
         * @brief
         *
         * @tparam ScalarIterable
         * @param stl_ranges
         */
        template <typename ScalarIterable>
        void read_from_scalar_array(ScalarIterable const &stl_ranges);

        template <typename ScalarIterable>
        void evaluate_general_derivative(ScalarIterable &stl_ranges);

        inline void collect_increment(bool sync) { sync_increment_ = sync; };

        void clear_increment() { calc::array_set_zero(increment_); };

        [[nodiscard]] size_t variable_number() const;

        // Friend functions
        template <CONCEPT_PARTICLES P, CONCEPT_INTERACTION F>
        friend std::ostream &operator<<(std::ostream &os, SimpleSystem<P, F> const &ps);

        template <CONCEPT_PARTICLES P, CONCEPT_INTERACTION F>
        friend std::istream &operator>>(std::istream &is, SimpleSystem<P, F> &ps);

        inline constexpr size_t time_offset() const { return 0; };
        inline constexpr size_t pos_offset() const { return 1; };
        inline constexpr size_t vel_offset() const { return this->number() * 3 + 1; };
        inline constexpr size_t auxi_vel_offset() const { return this->number() * 6 + 1; };

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

        template <typename Array>
        void sync_pos_increment(Array const &inc, Scalar step_size);

        template <typename Array>
        void sync_vel_increment(Array const &inc, Scalar step_size);

        template <typename Array>
        void sync_auxi_vel_increment(Array const &inc, Scalar step_size);

        void sync_time_increment(Scalar phy_time);

        // Private members
        // Particles ptcl_;

        force::InteractionData<Interactions, VectorArray> accels_;

        StateScalarArray increment_;

        std::conditional_t<Interactions::ext_vel_dep, StateVectorArray, Empty> aux_vel_;

        bool sync_increment_{false};
    };
}  // namespace space::particle_system

namespace space::particle_system {
    /*---------------------------------------------------------------------------*\
        Class SimpleSystem Implementation
    \*---------------------------------------------------------------------------*/
    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
    template <CONCEPT_PARTICLE_CONTAINER STL>
    SimpleSystem<Particles, Interactions>::SimpleSystem(Scalar time, const STL &particle_set)
        : Particles(time, particle_set), accels_(particle_set.size()), increment_(this->variable_number()) {
        if constexpr (Interactions::ext_vel_dep) {
            aux_vel_ = this->vel();
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
    template <typename ScalarIterable>
    void SimpleSystem<Particles, Interactions>::read_from_scalar_array(const ScalarIterable &stl_ranges) {
        if (stl_ranges.size() == this->variable_number()) {
            auto begin = stl_ranges.begin();
            this->time() = *(begin + time_offset());

            auto pos_begin = begin + pos_offset();
            auto pos_end = begin + vel_offset();
            auto vel_begin = begin + vel_offset();
            auto vel_end = begin + auxi_vel_offset();

            load_to_coords(pos_begin, pos_end, this->pos());
            load_to_coords(vel_begin, vel_end, this->vel());
            if constexpr (Interactions::ext_vel_dep) {
                auto aux_vel_begin = begin + auxi_vel_offset();
                auto aux_vel_end = stl_ranges.end();
                load_to_coords(aux_vel_begin, aux_vel_end, aux_vel_);
            }
        } else {
            spacehub_abort("Wrong input array size!");
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
    template <typename ScalarIterable>
    void SimpleSystem<Particles, Interactions>::write_to_scalar_array(ScalarIterable &stl_ranges) {
        stl_ranges.clear();
        stl_ranges.reserve(this->variable_number());
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
        stl_ranges.reserve(this->variable_number());
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
    void SimpleSystem<Particles, Interactions>::kick(Scalar step_size) {
        if constexpr (Interactions::ext_vel_dep) {
            Scalar half_step = 0.5 * step_size;

            eval_vel_indep_acc();

            kick_real_vel(half_step);
            kick_pseu_vel(step_size);
            kick_real_vel(half_step);
        } else {
            Interactions::eval_acc(*this, accels_.acc());
            calc::array_advance(this->vel(), accels_.acc(), step_size);
            sync_vel_increment(accels_.acc(), step_size);
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
    void SimpleSystem<Particles, Interactions>::drift(Scalar step_size) {
        this->time() += step_size;
        calc::array_advance(this->pos(), this->vel(), step_size);
        sync_time_increment(step_size);
        sync_pos_increment(this->vel(), step_size);
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
    template <typename GenVectorArray>
    void SimpleSystem<Particles, Interactions>::evaluate_acc(GenVectorArray &acceleration) const {
        Interactions::eval_acc(*this, acceleration);
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
    void SimpleSystem<Particles, Interactions>::kick_real_vel(Scalar step_size) {
        std::swap(aux_vel_, this->vel());
        Interactions::eval_extra_vel_dep_acc(*this, accels_.ext_vel_dep_acc());
        std::swap(aux_vel_, this->vel());
        calc::array_add(accels_.acc(), accels_.tot_vel_indep_acc(), accels_.ext_vel_dep_acc());
        calc::array_advance(this->vel(), accels_.acc(), step_size);
        sync_vel_increment(accels_.acc(), step_size);
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
    void SimpleSystem<Particles, Interactions>::kick_pseu_vel(Scalar step_size) {
        Interactions::eval_extra_vel_dep_acc(*this, accels_.ext_vel_dep_acc());
        calc::array_add(accels_.acc(), accels_.tot_vel_indep_acc(), accels_.ext_vel_dep_acc());
        calc::array_advance(aux_vel_, accels_.acc(), step_size);
        sync_auxi_vel_increment(accels_.acc(), step_size);
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
    template <typename Array>
    void SimpleSystem<Particles, Interactions>::sync_pos_increment(Array const &inc, Scalar step_size) {
        if (sync_increment_) {
            advance_scaled_coords_to(increment_.begin() + pos_offset(), inc, step_size);
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
    template <typename Array>
    void SimpleSystem<Particles, Interactions>::sync_vel_increment(Array const &inc, Scalar step_size) {
        if (sync_increment_) {
            advance_scaled_coords_to(increment_.begin() + vel_offset(), inc, step_size);
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
    template <typename Array>
    void SimpleSystem<Particles, Interactions>::sync_auxi_vel_increment(Array const &inc, Scalar step_size) {
        if (sync_increment_) {
            advance_scaled_coords_to(increment_.begin() + auxi_vel_offset(), inc, step_size);
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
    void SimpleSystem<Particles, Interactions>::sync_time_increment(Scalar phy_time) {
        if (sync_increment_) {
            increment_[time_offset()] += phy_time;
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
    std::ostream &operator<<(std::ostream &os, SimpleSystem<Particles, Interactions> const &ps) {
        os << static_cast<Particles>(ps);
        return os;
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions>
    std::istream &operator>>(std::istream &is, SimpleSystem<Particles, Interactions> &ps) {
        is >> static_cast<Particles>(ps);
        return is;
    }

}  // namespace space::particle_system
