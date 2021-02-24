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
 * @file archain.hpp
 *
 * Header file.
 */
#pragma once

#include "../type-class.hpp"
#include "chain.hpp"
#include "regu-system.hpp"
namespace space::particle_system {

    /*---------------------------------------------------------------------------*\
        Class ARchainSystem Declaration
    \*---------------------------------------------------------------------------*/
    /**
     *
     * @tparam Particles
     * @tparam Interactions
     * @tparam RegType
     */
    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType = ReguType::LogH>
    class ARchainSystem : public Particles {
       public:
        // Type members
        SPACEHUB_USING_TYPE_SYSTEM_OF(Particles);

        using Particle = typename Particles::Particle;

        using Interaction = Interactions;

        // Static public members
        static constexpr bool ext_vel_dep{Interactions::ext_vel_dep};

        static constexpr bool ext_vel_indep{Interactions::ext_vel_indep};

        static constexpr ReguType regu_type{RegType};

        // Constructors
        SPACEHUB_MAKE_CONSTRUCTORS(ARchainSystem, delete, default, default, default, default);

        ARchainSystem(Scalar time, concepts::ParticleContainer auto const &particle_set);

        // Public methods
        SPACEHUB_ARRAY_ACCESSOR(StateVectorArray, chain_pos, chain_pos_);

        SPACEHUB_ARRAY_ACCESSOR(StateVectorArray, chain_vel, chain_vel_);

        SPACEHUB_ARRAY_ACCESSOR(IdxArray, index, index_);

        SPACEHUB_STD_ACCESSOR(StateScalar, omega, regu_.omega());

        SPACEHUB_STD_ACCESSOR(StateScalar, bindE, regu_.bindE());

        SPACEHUB_ARRAY_ACCESSOR(StateScalarArray, increment, increment_);

        Scalar step_scale() const { return regu_.regu_function(*this); };

        void evaluate_acc(concepts::Vec3Container auto &acceleration) const;

        void drift(Scalar step_size);

        void kick(Scalar step_size);

        void pre_iter_process();

        void post_iter_process();

        void write_to_scalar_array(concepts::ScalarContainer auto &container);

        void read_from_scalar_array(concepts::ScalarContainer auto const &container);

        void clear_increment() { calc::array_set_zero(increment_); };

        size_t variable_number() const;

        void evaluate_general_derivative(concepts::ScalarContainer auto &container);

        inline void collect_increment(bool sync) { sync_increment_ = sync; };

        inline constexpr size_t time_offset() const { return 0; };

        inline constexpr size_t pos_offset() const { return 1; };

        inline constexpr size_t vel_offset() const { return this->number() * 3 + 1; };

        inline constexpr size_t auxi_vel_offset() const { return this->number() * 6 + 1; };

        inline constexpr size_t omega_offset() const {
            return this->number() * 3 * (2 + static_cast<size_t>(Interactions::ext_vel_dep)) + 1;
        };

        inline constexpr size_t bindE_offset() const {
            return this->number() * 3 * (2 + static_cast<size_t>(Interactions::ext_vel_dep)) + 2;
        };

        // Friend functions
        template <CONCEPT_PARTICLES P, CONCEPT_INTERACTION F, ReguType R>
        friend std::ostream &operator<<(std::ostream &os, ARchainSystem<P, F, R> const &ps);

        template <CONCEPT_PARTICLES P, CONCEPT_INTERACTION F, ReguType R>
        friend std::istream &operator>>(std::istream &is, ARchainSystem<P, F, R> &ps);

       private:
        // Private methods
        void chain_advance(concepts::Vec3Container auto &var, concepts::Vec3Container auto &chain_var,
                           concepts::Vec3Container auto const &chain_increment, Scalar phy_time);

        void eval_vel_indep_acc();

        void advance_omega(StateVectorArray const &velocity, VectorArray const &d_omega_dr, Scalar phy_time);

        void advance_bindE(StateVectorArray const &velocity, VectorArray const &d_bindE_dr, Scalar phy_time);

        void kick_pseu_vel(Scalar phy_time);

        void kick_real_vel(Scalar phy_time);

        void sync_pos_increment(concepts::Vec3Container auto const &inc, Scalar step_size);

        void sync_vel_increment(concepts::Vec3Container auto const &inc, Scalar step_size);

        void sync_auxi_vel_increment(concepts::Vec3Container auto const &inc, Scalar step_size);

        void sync_time_increment(Scalar phy_time);

        void sync_omega_increment(Scalar domega);

        void sync_bindE_increment(Scalar dbindE);

        // Private members
        force::InteractionData<Interactions, VectorArray> accels_;

        Regularization<TypeSet, RegType> regu_;

        StateVectorArray chain_pos_;

        StateVectorArray chain_vel_;

        VectorArray chain_acc_;

        IdxArray index_;

        IdxArray new_index_;

        StateScalarArray increment_;

        std::conditional_t<Interactions::ext_vel_dep, StateVectorArray, Empty> aux_vel_;

        std::conditional_t<Interactions::ext_vel_dep, StateVectorArray, Empty> chain_aux_vel_;

        bool sync_increment_{false};

        CREATE_MEMBER_CHECK(err);
    };

#define CLASS_ARchainSystem(...)                                                               \
    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType> \
    __VA_ARGS__ ARchainSystem<Particles, Interactions, RegType>

    /*---------------------------------------------------------------------------*\
        Class ARchainSystem Implementation
    \*---------------------------------------------------------------------------*/
    CLASS_ARchainSystem()::ARchainSystem(Scalar time, concepts::ParticleContainer auto const &particle_set)
        : Particles(time, particle_set),
          accels_(particle_set.size()),
          regu_(static_cast<Particles>(*this)),  // chain_pos that might be invoked by regu is not initialized yet.
          chain_pos_(particle_set.size()),
          chain_vel_(particle_set.size()),
          chain_acc_(particle_set.size()),
          index_(particle_set.size()),
          new_index_(particle_set.size()),
          increment_(this->variable_number()) {
        Chain::calc_chain_index(this->pos(), index_);
        Chain::calc_chain(this->pos(), chain_pos(), index());
        Chain::calc_chain(this->vel(), chain_vel(), index());
        if constexpr (Interactions::ext_vel_dep) {
            aux_vel_ = this->vel();
            chain_aux_vel_ = chain_vel_;
        }
        regu_ = std::move(
            Regularization<TypeSet, RegType>{*this});  // re-construct the regularization with chain coordinates
    }

    CLASS_ARchainSystem(void)::evaluate_acc(concepts::Vec3Container auto &acceleration) const {
        Interactions::eval_acc(*this, acceleration);
    }

    CLASS_ARchainSystem(void)::drift(Scalar step_size) {
        Scalar phy_time = regu_.eval_pos_phy_time(*this, step_size);
        chain_advance(this->pos(), chain_pos(), chain_vel(), phy_time);
        this->time() += phy_time;
        sync_time_increment(phy_time);
        sync_pos_increment(chain_vel(), phy_time);
    }

    CLASS_ARchainSystem(void)::kick(Scalar step_size) {
        Scalar phy_time = regu_.eval_vel_phy_time(*this, step_size);
        Scalar half_time = 0.5 * phy_time;

        eval_vel_indep_acc();

        if constexpr (Interactions::ext_vel_dep) {
            kick_real_vel(half_time);
            kick_pseu_vel(phy_time);
            kick_real_vel(half_time);
        } else {
            Chain::calc_chain(accels_.tot_vel_indep_acc(), chain_acc_, index());

            advance_omega(this->vel(), accels_.newtonian_acc(), half_time);
            if constexpr (Interactions::ext_vel_indep) {
                advance_bindE(this->vel(), accels_.ext_vel_indep_acc(), half_time);
            }
            chain_advance(this->vel(), chain_vel(), chain_acc_, phy_time);
            sync_vel_increment(chain_acc_, phy_time);
            if constexpr (Interactions::ext_vel_indep) {
                advance_bindE(this->vel(), accels_.ext_vel_indep_acc(), half_time);
            }
            advance_omega(this->vel(), accels_.newtonian_acc(), half_time);
        }
    }

    CLASS_ARchainSystem(void)::pre_iter_process() {
        if constexpr (Interactions::ext_vel_dep) {
            aux_vel_ = this->vel();
            chain_aux_vel_ = chain_vel_;
        }
    }

    CLASS_ARchainSystem(void)::post_iter_process() {
        Chain::calc_chain_index(this->pos(), new_index_);
        if (new_index_ != index_) {
            Chain::update_chain(chain_pos_, this->pos(), index_, new_index_);
            Chain::calc_cartesian(this->mass(), chain_pos_, this->pos(), new_index_);
            Chain::update_chain(chain_vel_, this->vel(), index_, new_index_);
            Chain::calc_cartesian(this->mass(), chain_vel_, this->vel(), new_index_);
            index_ = new_index_;
        }
    }

    CLASS_ARchainSystem(void)::write_to_scalar_array(concepts::ScalarContainer auto &container) {
        container.clear();
        container.reserve(this->variable_number());
        container.emplace_back(this->time());

        add_coords_to(container, chain_pos_);
        add_coords_to(container, chain_vel_);

        /*container.reserve(this->number() * 12 + 3);
        container.emplace_back(this->time());
        container.emplace_back(omega());
        container.emplace_back(bindE());
        add_coords_to(container, chain_pos_);
        add_coords_to(container, chain_vel_);
        add_coords_to(container, this->pos());
        add_coords_to(container, this->vel());*/
        if constexpr (Interactions::ext_vel_dep) {
            add_coords_to(container, chain_aux_vel_);
        }

        container.emplace_back(omega());
        container.emplace_back(bindE());
    }

    CLASS_ARchainSystem(void)::evaluate_general_derivative(concepts::ScalarContainer auto &container) {
        container.clear();
        container.reserve(this->variable_number());

        Scalar pos_regu = regu_.eval_pos_phy_time(*this, 1);
        Scalar vel_regu = regu_.eval_vel_phy_time(*this, 1);

        container.emplace_back(pos_regu);

        Interactions::eval_newtonian_acc(*this, accels_.newtonian_acc());
        if constexpr (Interactions::ext_vel_indep || Interactions::ext_vel_dep) {
            Interactions::eval_extra_acc(*this, accels_.acc());
            calc::array_add(accels_.acc(), accels_.acc(), accels_.newtonian_acc());
        }

        add_scaled_coords_to(container, chain_vel_, pos_regu);

        if constexpr (Interactions::ext_vel_indep || Interactions::ext_vel_dep) {
            Chain::calc_chain(accels_.acc(), chain_acc_, index());
            add_scaled_coords_to(container, chain_acc_, vel_regu);
        } else {
            Chain::calc_chain(accels_.newtonian_acc(), chain_acc_, index());
            add_scaled_coords_to(container, chain_acc_, vel_regu);
        }
        if constexpr (Interactions::ext_vel_dep) {
            add_scaled_coords_to(container, chain_acc_, vel_regu);
        }

        if constexpr (regu_type == ReguType::TTL) {
            Scalar d_omega_dh =
                calc::coord_contract_to_scalar(this->mass(), this->vel(), accels_.newtonian_acc()) * vel_regu;
            container.emplace_back(d_omega_dh);
        } else {
            container.emplace_back(0);
        }

        if constexpr (Interactions::ext_vel_indep || Interactions::ext_vel_dep) {
            if constexpr (regu_type == ReguType::LogH) {
                Scalar d_bindE_dh =
                    -calc::coord_contract_to_scalar(this->mass(), this->vel(), accels_.acc()) * vel_regu;
                container.emplace_back(d_bindE_dh);
            } else {
                container.emplace_back(0);
            }

        } else {
            container.emplace_back(0);
        }
    }

    CLASS_ARchainSystem(void)::sync_pos_increment(concepts::Vec3Container auto const &inc, Scalar step_size) {
        if (sync_increment_) {
            size_t offset = pos_offset();
            for (size_t i = 0; i < inc.size(); ++i) {
                increment_[3 * i + offset] += inc[i].x * step_size;
                increment_[3 * i + offset + 1] += inc[i].y * step_size;
                increment_[3 * i + offset + 2] += inc[i].z * step_size;
            }
        }
    }

    CLASS_ARchainSystem(void)::sync_vel_increment(concepts::Vec3Container auto const &inc, Scalar step_size) {
        if (sync_increment_) {
            size_t offset = vel_offset();
            for (size_t i = 0; i < inc.size(); ++i) {
                increment_[3 * i + offset] += inc[i].x * step_size;
                increment_[3 * i + offset + 1] += inc[i].y * step_size;
                increment_[3 * i + offset + 2] += inc[i].z * step_size;
            }
        }
    }

    CLASS_ARchainSystem(void)::sync_auxi_vel_increment(concepts::Vec3Container auto const &inc, Scalar step_size) {
        if (sync_increment_) {
            size_t offset = auxi_vel_offset();
            for (size_t i = 0; i < inc.size(); ++i) {
                increment_[3 * i + offset] += inc[i].x * step_size;
                increment_[3 * i + offset + 1] += inc[i].y * step_size;
                increment_[3 * i + offset + 2] += inc[i].z * step_size;
            }
        }
    }

    CLASS_ARchainSystem(void)::sync_time_increment(Scalar phy_time) {
        if (sync_increment_) {
            increment_[time_offset()] += phy_time;
        }
    }

    CLASS_ARchainSystem(void)::sync_omega_increment(Scalar domega) {
        if (sync_increment_) {
            increment_[omega_offset()] += domega;
        }
    }

    CLASS_ARchainSystem(void)::sync_bindE_increment(Scalar dbindE) {
        if (sync_increment_) {
            increment_[bindE_offset()] += dbindE;
        }
    }

    CLASS_ARchainSystem(void)::read_from_scalar_array(const concepts::ScalarContainer auto &container) {
        if (container.size() == this->variable_number()) {
            auto begin = container.begin();
            this->time() = *(begin + time_offset());

            auto pos_begin = begin + pos_offset();
            auto pos_end = begin + vel_offset();
            auto vel_begin = begin + vel_offset();
            auto vel_end = begin + auxi_vel_offset();
            load_to_coords(pos_begin, pos_end, chain_pos_);
            load_to_coords(vel_begin, vel_end, chain_vel_);

            /*auto c_pos_begin = vel_end;
            auto c_pos_end = c_pos_begin + len;
            auto c_vel_begin = c_pos_end;
            auto c_vel_end = c_vel_begin + len;

            load_to_coords(c_pos_begin, c_pos_end, this->pos());
            load_to_coords(c_vel_begin, c_vel_end, this->vel());*/

            Chain::calc_cartesian(this->mass(), chain_pos_, this->pos(), index());
            Chain::calc_cartesian(this->mass(), chain_vel_, this->vel(), index());

            if constexpr (Interactions::ext_vel_dep) {
                auto aux_vel_begin = begin + auxi_vel_offset();
                auto aux_vel_end = container.end();
                load_to_coords(aux_vel_begin, aux_vel_end, chain_aux_vel_);
                Chain::calc_cartesian(this->mass(), chain_aux_vel_, aux_vel_, index_);
            }

            omega() = *(begin + omega_offset());
            bindE() = *(begin + bindE_offset());
        } else {
            spacehub_abort("Wrong input array size!");
        }
    }

    CLASS_ARchainSystem(size_t)::variable_number() const {
        return this->number() * 3 * (2 + static_cast<size_t>(Interactions::ext_vel_dep)) + 3;
    }

    CLASS_ARchainSystem(void)::chain_advance(concepts::Vec3Container auto &var, concepts::Vec3Container auto &chain_var,
                                             concepts::Vec3Container auto const &chain_increment, Scalar phy_time) {
        calc::array_advance(chain_var, chain_increment, phy_time);
        Chain::calc_cartesian(this->mass(), chain_var, var, index());
    }

    CLASS_ARchainSystem(void)::eval_vel_indep_acc() {
        Interactions::eval_newtonian_acc(*this, accels_.newtonian_acc());

        if constexpr (Interactions::ext_vel_indep) {
            Interactions::eval_extra_vel_indep_acc(*this, accels_.ext_vel_indep_acc());
            calc::array_add(accels_.tot_vel_indep_acc(), accels_.ext_vel_indep_acc(), accels_.newtonian_acc());
        } else {
            accels_.tot_vel_indep_acc() = accels_.newtonian_acc();
        }
    }

    CLASS_ARchainSystem(void)::advance_omega(StateVectorArray const &velocity, VectorArray const &d_omega_dr,
                                             Scalar phy_time) {
        if constexpr (regu_type == ReguType::TTL) {
            Scalar d_omega = calc::coord_contract_to_scalar(this->mass(), velocity, d_omega_dr) * phy_time;
            regu_.omega() += d_omega;
            sync_omega_increment(d_omega);
        }
    }

    CLASS_ARchainSystem(void)::advance_bindE(StateVectorArray const &velocity, VectorArray const &d_bindE_dr,
                                             Scalar phy_time) {
        if constexpr ((Interactions::ext_vel_indep || Interactions::ext_vel_dep) && regu_type == ReguType::LogH) {
            Scalar d_bindE = -calc::coord_contract_to_scalar(this->mass(), velocity, d_bindE_dr) * phy_time;
            regu_.bindE() += d_bindE;
            sync_bindE_increment(d_bindE);
        }
    }

    CLASS_ARchainSystem(void)::kick_pseu_vel(Scalar phy_time) {
        Interactions::eval_extra_vel_dep_acc(*this, accels_.ext_vel_dep_acc());
        calc::array_add(accels_.acc(), accels_.tot_vel_indep_acc(), accels_.ext_vel_dep_acc());
        Chain::calc_chain(accels_.acc(), chain_acc_, index());
        chain_advance(aux_vel_, chain_aux_vel_, chain_acc_, phy_time);
        sync_auxi_vel_increment(chain_acc_, phy_time);
        advance_omega(this->vel(), accels_.newtonian_acc(), phy_time);

        if constexpr (Interactions::ext_vel_indep) {
            calc::array_add(accels_.acc(), accels_.ext_vel_indep_acc(), accels_.ext_vel_dep_acc());
            advance_bindE(this->vel(), accels_.acc(), phy_time);
        } else {
            advance_bindE(this->vel(), accels_.ext_vel_dep_acc(), phy_time);
        }
    }

    CLASS_ARchainSystem(void)::kick_real_vel(Scalar phy_time) {
        std::swap(aux_vel_, this->vel());
        std::swap(chain_aux_vel_, chain_vel());
        Interactions::eval_extra_vel_dep_acc(*this, accels_.ext_vel_dep_acc());
        std::swap(aux_vel_, this->vel());
        std::swap(chain_aux_vel_, chain_vel());
        calc::array_add(accels_.acc(), accels_.tot_vel_indep_acc(), accels_.ext_vel_dep_acc());

        Chain::calc_chain(accels_.acc(), chain_acc_, index());
        chain_advance(this->vel(), chain_vel(), chain_acc_, phy_time);
        sync_vel_increment(chain_acc_, phy_time);
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    std::ostream &operator<<(std::ostream &os, ARchainSystem<Particles, Interactions, RegType> const &ps) {
        os << static_cast<Particles>(ps);
        return os;
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    std::istream &operator>>(std::istream &is, ARchainSystem<Particles, Interactions, RegType> &ps) {
        is >> static_cast<Particles>(ps);
        return is;
    }
}  // namespace space::particle_system
