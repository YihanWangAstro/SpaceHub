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

        SimpleSystem(Scalar time, concepts::ParticleContainer auto const &particle_set);

        // Public methods
        SPACEHUB_ARRAY_ACCESSOR(StateScalarArray, increment, increment_);

        Scalar step_scale() const { return 1.0; };

        void evaluate_acc(concepts::Vec3Container auto &acceleration) const;

        void drift(Scalar step_size);

        void kick(Scalar step_size);

        void pre_iter_process();

        void post_iter_process(){};

        void write_to_scalar_array(concepts::ScalarContainer auto &container);

        void read_from_scalar_array(concepts::ScalarContainer auto const &container);

        void evaluate_general_derivative(concepts::ScalarContainer auto &container);

        inline void collect_increment(bool sync) { sync_increment_ = sync; };

        void clear_increment() { calc::array_set_zero(increment_); };

        size_t variable_number() const;

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
        void eval_vel_indep_acc();

        void kick_pseu_vel(Scalar step_size);

        void kick_real_vel(Scalar step_size);

        void sync_pos_increment(concepts::Vec3Container auto const &inc, Scalar step_size);

        void sync_vel_increment(concepts::Vec3Container auto const &inc, Scalar step_size);

        void sync_auxi_vel_increment(concepts::Vec3Container auto const &inc, Scalar step_size);

        void sync_time_increment(Scalar phy_time);

        // Private members
        force::InteractionData<Interactions, VectorArray> accels_;

        StateScalarArray increment_;

        std::conditional_t<Interactions::ext_vel_dep, StateVectorArray, Empty> aux_vel_;

        bool sync_increment_{false};
    };
}  // namespace space::particle_system

namespace space::particle_system {

#define CLASS_SimpleSystem(...)                                              \
    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions> \
    __VA_ARGS__ SimpleSystem<Particles, Interactions>

    /*---------------------------------------------------------------------------*\
        Class SimpleSystem Implementation
    \*---------------------------------------------------------------------------*/
    CLASS_SimpleSystem()::SimpleSystem(Scalar time, concepts::ParticleContainer auto const &particle_set)
        : Particles(time, particle_set), accels_(particle_set.size()), increment_(this->variable_number()) {
        if constexpr (Interactions::ext_vel_dep) {
            aux_vel_ = this->vel();
        }
    }

    CLASS_SimpleSystem(void)::read_from_scalar_array(concepts::ScalarContainer auto const &container) {
        if (container.size() == this->variable_number()) {
            auto begin = container.begin();
            this->time() = *(begin + time_offset());

            auto pos_begin = begin + pos_offset();
            auto pos_end = begin + vel_offset();
            auto vel_begin = begin + vel_offset();
            auto vel_end = begin + auxi_vel_offset();

            load_to_coords(pos_begin, pos_end, this->pos());
            load_to_coords(vel_begin, vel_end, this->vel());
            if constexpr (Interactions::ext_vel_dep) {
                auto aux_vel_begin = begin + auxi_vel_offset();
                auto aux_vel_end = container.end();
                load_to_coords(aux_vel_begin, aux_vel_end, aux_vel_);
            }
        } else {
            spacehub_abort("Wrong input array size!");
        }
    }

    CLASS_SimpleSystem(void)::write_to_scalar_array(concepts::ScalarContainer auto &container) {
        container.clear();
        container.reserve(this->variable_number());
        container.emplace_back(this->time());
        add_coords_to(container, this->pos());
        add_coords_to(container, this->vel());
        if constexpr (Interactions::ext_vel_dep) {
            add_coords_to(container, aux_vel_);
        }
    }

    CLASS_SimpleSystem(void)::evaluate_general_derivative(concepts::ScalarContainer auto &container) {
        container.clear();
        container.reserve(this->variable_number());
        container.emplace_back(1);              // dt/dh
        add_coords_to(container, this->vel());  // dp/dt
        Interactions::eval_acc(*this, this->accels_.acc());
        add_coords_to(container, this->accels_.acc());  // dv/dt
        if constexpr (Interactions::ext_vel_dep) {
            add_coords_to(container, this->accels_.acc());  // dw/dt
        }
    }

    CLASS_SimpleSystem(size_t)::variable_number() const {
        return this->number() * 3 * (2 + static_cast<size_t>(Interactions::ext_vel_dep)) + 1;
    }

    CLASS_SimpleSystem(void)::pre_iter_process() {
        if constexpr (Interactions::ext_vel_dep) {
            aux_vel_ = this->vel();
        }
    }

    CLASS_SimpleSystem(void)::kick(Scalar step_size) {
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

    CLASS_SimpleSystem(void)::drift(Scalar step_size) {
        this->time() += step_size;
        calc::array_advance(this->pos(), this->vel(), step_size);
        sync_time_increment(step_size);
        sync_pos_increment(this->vel(), step_size);
    }

    CLASS_SimpleSystem(void)::evaluate_acc(concepts::Vec3Container auto &acceleration) const {
        Interactions::eval_acc(*this, acceleration);
    }

    CLASS_SimpleSystem(void)::kick_real_vel(Scalar step_size) {
        std::swap(aux_vel_, this->vel());
        Interactions::eval_extra_vel_dep_acc(*this, accels_.ext_vel_dep_acc());
        std::swap(aux_vel_, this->vel());
        calc::array_add(accels_.acc(), accels_.tot_vel_indep_acc(), accels_.ext_vel_dep_acc());
        calc::array_advance(this->vel(), accels_.acc(), step_size);
        sync_vel_increment(accels_.acc(), step_size);
    }

    CLASS_SimpleSystem(void)::kick_pseu_vel(Scalar step_size) {
        Interactions::eval_extra_vel_dep_acc(*this, accels_.ext_vel_dep_acc());
        calc::array_add(accels_.acc(), accels_.tot_vel_indep_acc(), accels_.ext_vel_dep_acc());
        calc::array_advance(aux_vel_, accels_.acc(), step_size);
        sync_auxi_vel_increment(accels_.acc(), step_size);
    }

    CLASS_SimpleSystem(void)::eval_vel_indep_acc() {
        Interactions::eval_newtonian_acc(*this, accels_.tot_vel_indep_acc());
        if constexpr (Interactions::ext_vel_indep) {
            Interactions::eval_extra_vel_indep_acc(*this, accels_.ext_vel_indep_acc());
            calc::array_add(accels_.tot_vel_indep_acc(), accels_.tot_vel_indep_acc(), accels_.ext_vel_indep_acc());
        }
    }

    CLASS_SimpleSystem(void)::sync_pos_increment(concepts::Vec3Container auto const &inc, Scalar step_size) {
        if (sync_increment_) {
            size_t offset = pos_offset();
            for (size_t i = 0; i < inc.size(); ++i) {
                increment_[3 * i + offset] += inc[i].x * step_size;
                increment_[3 * i + offset + 1] += inc[i].y * step_size;
                increment_[3 * i + offset + 2] += inc[i].z * step_size;
            }
        }
    }

    CLASS_SimpleSystem(void)::sync_vel_increment(concepts::Vec3Container auto const &inc, Scalar step_size) {
        if (sync_increment_) {
            size_t offset = vel_offset();
            for (size_t i = 0; i < inc.size(); ++i) {
                increment_[3 * i + offset] += inc[i].x * step_size;
                increment_[3 * i + offset + 1] += inc[i].y * step_size;
                increment_[3 * i + offset + 2] += inc[i].z * step_size;
            }
        }
    }

    CLASS_SimpleSystem(void)::sync_auxi_vel_increment(concepts::Vec3Container auto const &inc, Scalar step_size) {
        if (sync_increment_) {
            size_t offset = auxi_vel_offset();
            for (size_t i = 0; i < inc.size(); ++i) {
                increment_[3 * i + offset] += inc[i].x * step_size;
                increment_[3 * i + offset + 1] += inc[i].y * step_size;
                increment_[3 * i + offset + 2] += inc[i].z * step_size;
            }
        }
    }

    CLASS_SimpleSystem(void)::sync_time_increment(Scalar phy_time) {
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
