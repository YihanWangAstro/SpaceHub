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
    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    class ARchainSystem : public Particles {
       public:
        // Type members
        SPACEHUB_USING_TYPE_SYSTEM_OF(Particles);

        using Particle = typename Particles::Particle;

        using Interaction = Interactions;

        // Static public members
        static constexpr ReguType regu_type{RegType};

        // Constructors
        SPACEHUB_MAKE_CONSTRUCTORS(ARchainSystem, delete, default, default, default, default);

        template <CONCEPT_PARTICLE_CONTAINER STL>
        ARchainSystem(Scalar time, STL const &particle_set);

        // Public methods
        SPACEHUB_STD_ACCESSOR(auto, chain_pos, chain_pos_);

        SPACEHUB_STD_ACCESSOR(auto, chain_vel, chain_vel_);

        SPACEHUB_STD_ACCESSOR(auto, index, index_);

        SPACEHUB_STD_ACCESSOR(auto, omega, regu_.omega());

        SPACEHUB_STD_ACCESSOR(auto, bindE, regu_.bindE());

        Scalar regu_function() const { return regu_.regu_function(*this); };

        void advance_time(Scalar step_size);

        template <typename GenVectorArray>
        void advance_pos(Scalar step_size, GenVectorArray const &velocity);

        template <typename GenVectorArray>
        void advance_vel(Scalar step_size, GenVectorArray const &acceleration);

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

        [[nodiscard]] size_t variable_number() const;

        /**
         * @brief
         *
         * @tparam STL
         * @param stl_ranges
         */
        template <typename ScalarIterable>
        void evaluate_general_derivative(ScalarIterable &stl_ranges);

        // Friend functions
        template <CONCEPT_PARTICLES P, CONCEPT_INTERACTION F, ReguType R>
        friend std::ostream &operator<<(std::ostream &os, ARchainSystem<P, F, R> const &ps);

        template <CONCEPT_PARTICLES P, CONCEPT_INTERACTION F, ReguType R>
        friend std::istream &operator>>(std::istream &is, ARchainSystem<P, F, R> &ps);

       private:
        // Private methods
        template <typename Array1, typename Array2, typename Array3>
        void chain_advance(Array1 &var, Array2 &chain_var, Array3 const &chain_increment, Scalar phy_time);

        void eval_vel_indep_acc();

        void advance_omega(AdVectorArray const &velocity, VectorArray const &d_omega_dr, Scalar phy_time);

        void advance_bindE(AdVectorArray const &velocity, VectorArray const &d_bindE_dr, Scalar phy_time);

        void kick_pseu_vel(Scalar phy_time);

        void kick_real_vel(Scalar phy_time);

        // Private members
        interactions::InteractionData<Interactions, VectorArray> accels_;

        Regularization<TypeSet, RegType> regu_;

        AdVectorArray chain_pos_;

        AdVectorArray chain_vel_;

        VectorArray chain_acc_;

        IdxArray index_;

        IdxArray new_index_;

        std::conditional_t<Interactions::ext_vel_dep, AdVectorArray, Empty> aux_vel_;

        std::conditional_t<Interactions::ext_vel_dep, AdVectorArray, Empty> chain_aux_vel_;
    };

    /*---------------------------------------------------------------------------*\
        Class ARchainSystem Implementation
    \*---------------------------------------------------------------------------*/
    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    template <CONCEPT_PARTICLE_CONTAINER STL>
    ARchainSystem<Particles, Interactions, RegType>::ARchainSystem(Scalar time, const STL &particle_set)
        : Particles(time, particle_set),
          accels_(particle_set.size()),
          regu_(static_cast<Particles>(*this)),  // chain_pos that might be invoked by regu is not initialized yet.
          chain_pos_(particle_set.size()),
          chain_vel_(particle_set.size()),
          chain_acc_(particle_set.size()),
          index_(particle_set.size()),
          new_index_(particle_set.size()) {
        Chain::calc_chain_index(this->pos(), index_);
        Chain::calc_chain(this->pos(), chain_pos(), index());
        Chain::calc_chain(this->vel(), chain_vel(), index());
        if constexpr (Interactions::ext_vel_dep) {
            aux_vel_ = this->vel();
            chain_aux_vel_ = chain_vel_;
        }
        regu_ = std::move(
            Regularization<TypeSet, RegType>{*this});  // re-construct the regularization with chain coordinates.
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    void ARchainSystem<Particles, Interactions, RegType>::advance_time(Scalar step_size) {
        Scalar phy_time = regu_.eval_pos_phy_time(*this, step_size);
        this->time() += phy_time;
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    template <typename GenVectorArray>
    void ARchainSystem<Particles, Interactions, RegType>::advance_pos(Scalar step_size,
                                                                      GenVectorArray const &velocity) {
        Scalar phy_time = regu_.eval_pos_phy_time(*this, step_size);
        Chain::calc_chain(velocity, chain_acc_, index());  // borrow chain_acc_ as chain velocity buffer.
        chain_advance(this->pos(), chain_pos(), chain_acc_, phy_time);
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    template <typename GenVectorArray>
    void ARchainSystem<Particles, Interactions, RegType>::advance_vel(Scalar step_size,
                                                                      GenVectorArray const &acceleration) {
        Scalar phy_time = regu_.eval_vel_phy_time(*this, step_size);
        Scalar half_time = 0.5 * phy_time;
        Chain::calc_chain(acceleration, chain_acc_, index());

        if constexpr (regu_type == ReguType::TTL) {
            chain_advance(this->vel(), chain_vel(), chain_acc_, half_time);
            Interactions::eval_newtonian_acc(*this, accels_.newtonian_acc());
            advance_omega(this->vel(), accels_.newtonian_acc(), phy_time);
            chain_advance(this->vel(), chain_vel(), chain_acc_, half_time);
        } else if constexpr (regu_type == ReguType::LogH) {
            if constexpr (Interactions::ext_vel_indep || Interactions::ext_vel_dep) {
                chain_advance(this->vel(), chain_vel(), chain_acc_, half_time);
                Interactions::eval_extra_acc(*this, accels_.acc());
                advance_bindE(this->vel(), accels_.acc(), phy_time);
                chain_advance(this->vel(), chain_vel(), chain_acc_, half_time);
            } else {
                chain_advance(this->vel(), chain_vel(), chain_acc_, phy_time);
            }
        } else {
            chain_advance(this->vel(), chain_vel(), chain_acc_, phy_time);
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    template <typename GenVectorArray>
    void ARchainSystem<Particles, Interactions, RegType>::evaluate_acc(GenVectorArray &acceleration) const {
        Interactions::eval_acc(*this, acceleration);
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    void ARchainSystem<Particles, Interactions, RegType>::drift(Scalar step_size) {
        Scalar phy_time = regu_.eval_pos_phy_time(*this, step_size);
        chain_advance(this->pos(), chain_pos(), chain_vel(), phy_time);
        this->time() += phy_time;
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    void ARchainSystem<Particles, Interactions, RegType>::kick(Scalar step_size) {
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
            if constexpr (Interactions::ext_vel_indep) {
                advance_bindE(this->vel(), accels_.ext_vel_indep_acc(), half_time);
            }
            advance_omega(this->vel(), accels_.newtonian_acc(), half_time);
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    void ARchainSystem<Particles, Interactions, RegType>::pre_iter_process() {
        if constexpr (Interactions::ext_vel_dep) {
            aux_vel_ = this->vel();
            chain_aux_vel_ = chain_vel_;
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    void ARchainSystem<Particles, Interactions, RegType>::post_iter_process() {
        Chain::calc_chain_index(this->pos(), new_index_);
        if (new_index_ != index_) {
            Chain::update_chain(chain_pos_, this->pos(), index_, new_index_);
            Chain::calc_cartesian(this->mass(), chain_pos_, this->pos(), new_index_);
            Chain::update_chain(chain_vel_, this->vel(), index_, new_index_);
            Chain::calc_cartesian(this->mass(), chain_vel_, this->vel(), new_index_);
            index_ = new_index_;
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    template <typename ScalarIterable>
    void ARchainSystem<Particles, Interactions, RegType>::write_to_scalar_array(ScalarIterable &stl_ranges) {
        stl_ranges.clear();
        stl_ranges.reserve(this->number() * 3 * (2 + static_cast<size_t>(Interactions::ext_vel_dep)) + 3);
        stl_ranges.emplace_back(this->time());
        stl_ranges.emplace_back(omega());
        stl_ranges.emplace_back(bindE());
        add_coords_to(stl_ranges, chain_pos_);
        add_coords_to(stl_ranges, chain_vel_);

        /*stl_ranges.reserve(this->number() * 12 + 3);
        stl_ranges.emplace_back(this->time());
        stl_ranges.emplace_back(omega());
        stl_ranges.emplace_back(bindE());
        add_coords_to(stl_ranges, chain_pos_);
        add_coords_to(stl_ranges, chain_vel_);
        add_coords_to(stl_ranges, this->pos());
        add_coords_to(stl_ranges, this->vel());*/
        if constexpr (Interactions::ext_vel_dep) {
            add_coords_to(stl_ranges, chain_aux_vel_);
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    template <typename ScalarIterable>
    void ARchainSystem<Particles, Interactions, RegType>::evaluate_general_derivative(ScalarIterable &stl_ranges) {
        stl_ranges.clear();
        stl_ranges.reserve(this->number() * 3 * (2 + static_cast<size_t>(Interactions::ext_vel_dep)) + 3);

        Scalar pos_regu = regu_.eval_pos_phy_time(*this, 1);
        Scalar vel_regu = regu_.eval_vel_phy_time(*this, 1);

        stl_ranges.emplace_back(pos_regu);

        Interactions::eval_newtonian_acc(*this, accels_.newtonian_acc());

        if constexpr (regu_type == ReguType::TTL) {
            Scalar d_omega_dh =
                calc::coord_contract_to_scalar(this->mass(), this->vel(), accels_.newtonian_acc()) * vel_regu;
            stl_ranges.emplace_back(d_omega_dh);
        } else {
            stl_ranges.emplace_back(0);
        }

        if constexpr (Interactions::ext_vel_indep || Interactions::ext_vel_dep) {
            Interactions::eval_extra_acc(*this, accels_.acc());
            if constexpr (regu_type == ReguType::LogH) {
                Scalar d_bindE_dh =
                    -calc::coord_contract_to_scalar(this->mass(), this->vel(), accels_.acc()) * vel_regu;
                stl_ranges.emplace_back(d_bindE_dh);
            } else {
                stl_ranges.emplace_back(0);
            }
            calc::array_add(accels_.acc(), accels_.acc(), accels_.newtonian_acc());
        } else {
            stl_ranges.emplace_back(0);
        }

        add_scaled_coords_to(stl_ranges, chain_vel_, pos_regu);

        if constexpr (Interactions::ext_vel_indep || Interactions::ext_vel_dep) {
            Chain::calc_chain(accels_.acc(), chain_acc_, index());
            add_scaled_coords_to(stl_ranges, chain_acc_, vel_regu);
        } else {
            Chain::calc_chain(accels_.newtonian_acc(), chain_acc_, index());
            add_scaled_coords_to(stl_ranges, chain_acc_, vel_regu);
        }
        if constexpr (Interactions::ext_vel_dep) {
            add_scaled_coords_to(stl_ranges, chain_acc_, vel_regu);
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    template <typename ScalarIterable>
    void ARchainSystem<Particles, Interactions, RegType>::read_from_scalar_array(const ScalarIterable &stl_ranges) {
        if (stl_ranges.size() == this->variable_number()) {
            auto begin = stl_ranges.begin();
            this->time() = *begin;
            omega() = *(begin + 1);
            bindE() = *(begin + 2);

            size_t len = this->number() * 3;

            auto pos_begin = begin + 3;
            auto pos_end = pos_begin + len;
            auto vel_begin = pos_end;
            auto vel_end = vel_begin + len;
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
                auto aux_vel_begin = vel_end;
                auto aux_vel_end = aux_vel_begin + len;
                load_to_coords(aux_vel_begin, aux_vel_end, chain_aux_vel_);
                Chain::calc_cartesian(this->mass(), chain_aux_vel_, aux_vel_, index_);
            }
        } else {
            spacehub_abort("Wrong input array size!");
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    size_t ARchainSystem<Particles, Interactions, RegType>::variable_number() const {
        return this->number() * 3 * (2 + static_cast<size_t>(Interactions::ext_vel_dep)) + 3;
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

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    template <typename Array1, typename Array2, typename Array3>
    void ARchainSystem<Particles, Interactions, RegType>::chain_advance(Array1 &var, Array2 &chain_var,
                                                                        Array3 const &chain_increment,
                                                                        Scalar phy_time) {
        calc::array_advance(chain_var, chain_increment, phy_time);
        Chain::calc_cartesian(this->mass(), chain_var, var, index());
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    void ARchainSystem<Particles, Interactions, RegType>::eval_vel_indep_acc() {
        Interactions::eval_newtonian_acc(*this, accels_.newtonian_acc());

        if constexpr (Interactions::ext_vel_indep) {
            Interactions::eval_extra_vel_indep_acc(*this, accels_.ext_vel_indep_acc());
            calc::array_add(accels_.tot_vel_indep_acc(), accels_.ext_vel_indep_acc(), accels_.newtonian_acc());
        } else {
            accels_.tot_vel_indep_acc() = accels_.newtonian_acc();
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    void ARchainSystem<Particles, Interactions, RegType>::advance_omega(AdVectorArray const &velocity,
                                                                        VectorArray const &d_omega_dr,
                                                                        Scalar phy_time) {
        if constexpr (regu_type == ReguType::TTL) {
            Scalar d_omega = calc::coord_contract_to_scalar(this->mass(), velocity, d_omega_dr);
            regu_.omega() += d_omega * phy_time;
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    void ARchainSystem<Particles, Interactions, RegType>::advance_bindE(AdVectorArray const &velocity,
                                                                        VectorArray const &d_bindE_dr,
                                                                        Scalar phy_time) {
        if constexpr ((Interactions::ext_vel_indep || Interactions::ext_vel_dep) && regu_type == ReguType::LogH) {
            Scalar d_bindE = -calc::coord_contract_to_scalar(this->mass(), velocity, d_bindE_dr);
            regu_.bindE() += d_bindE * phy_time;
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    void ARchainSystem<Particles, Interactions, RegType>::kick_pseu_vel(Scalar phy_time) {
        Interactions::eval_extra_vel_dep_acc(*this, accels_.ext_vel_dep_acc());
        calc::array_add(accels_.acc(), accels_.tot_vel_indep_acc(), accels_.ext_vel_dep_acc());
        Chain::calc_chain(accels_.acc(), chain_acc_, index());
        chain_advance(aux_vel_, chain_aux_vel_, chain_acc_, phy_time);
        advance_omega(this->vel(), accels_.newtonian_acc(), phy_time);

        if constexpr (Interactions::ext_vel_indep) {
            calc::array_add(accels_.acc(), accels_.ext_vel_indep_acc(), accels_.ext_vel_dep_acc());
            advance_bindE(this->vel(), accels_.acc(), phy_time);
        } else {
            advance_bindE(this->vel(), accels_.ext_vel_dep_acc(), phy_time);
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    void ARchainSystem<Particles, Interactions, RegType>::kick_real_vel(Scalar phy_time) {
        std::swap(aux_vel_, this->vel());
        std::swap(chain_aux_vel_, chain_vel());
        Interactions::eval_extra_vel_dep_acc(*this, accels_.ext_vel_dep_acc());
        std::swap(aux_vel_, this->vel());
        std::swap(chain_aux_vel_, chain_vel());
        calc::array_add(accels_.acc(), accels_.tot_vel_indep_acc(), accels_.ext_vel_dep_acc());

        Chain::calc_chain(accels_.acc(), chain_acc_, index());
        chain_advance(this->vel(), chain_vel(), chain_acc_, phy_time);
    }
}  // namespace space::particle_system
