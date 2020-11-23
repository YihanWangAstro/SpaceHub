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
    class ARchainSystem {
       public:
        // Type members
        SPACEHUB_USING_TYPE_SYSTEM_OF(Particles);

        using Particle = typename Particles::Particle;

        // Static public members
        static constexpr ReguType regu_type{RegType};

        // Constructors
        SPACEHUB_MAKE_CONSTRUCTORS(ARchainSystem, delete, default, default, default, default);

        template <CONCEPT_PARTICLE_CONTAINER STL>
        ARchainSystem(Scalar time, STL const &particle_set);

        // Public methods
        SPACEHUB_READ_ACCESSOR(Particles, particles, ptcl_);

        SPACEHUB_READ_ACCESSOR(Scalar, time, ptcl_.time());

        SPACEHUB_ARRAY_READ_ACCESSOR(IdxArray, idn, ptcl_.idn());

        SPACEHUB_ARRAY_READ_ACCESSOR(ScalarArray, mass, ptcl_.mass());

        SPACEHUB_ARRAY_READ_ACCESSOR(VectorArray, pos, ptcl_.pos());

        SPACEHUB_ARRAY_READ_ACCESSOR(VectorArray, vel, ptcl_.vel());

        SPACEHUB_STD_ACCESSOR(auto, chain_pos, chain_pos_);

        SPACEHUB_STD_ACCESSOR(auto, chain_vel, chain_vel_);

        SPACEHUB_STD_ACCESSOR(auto, index, index_);

        SPACEHUB_STD_ACCESSOR(auto, omega, regu_.omega());

        SPACEHUB_STD_ACCESSOR(auto, bindE, regu_.bindE());

        size_t number() const { return ptcl_.number(); };

        void advance_time(Scalar step_size);

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

        std::string column_names() const;

        // Friend functions
        template <CONCEPT_PARTICLES P, CONCEPT_INTERACTION F, ReguType R>
        friend std::ostream &operator<<(std::ostream &os, ARchainSystem<P, F, R> const &ps);

        template <CONCEPT_PARTICLES P, CONCEPT_INTERACTION F, ReguType R>
        friend std::istream &operator>>(std::istream &is, ARchainSystem<P, F, R> &ps);

       private:
        // Private methods
        void chain_advance(VectorArray &var, VectorArray &chain_var, VectorArray const &chain_increment,
                           Scalar phy_time);

        void eval_vel_indep_acc();

        void advance_omega(VectorArray const &velocity, VectorArray const &d_omega_dr, Scalar phy_time);

        void advance_bindE(VectorArray const &velocity, VectorArray const &d_bindE_dr, Scalar phy_time);

        void kick_pseu_vel(Scalar phy_time);

        void kick_real_vel(Scalar phy_time);

        // Private members
        Particles ptcl_;

        interactions::InteractionData<Interactions, VectorArray> accels_;

        Regularization<Scalar, RegType> regu_;

        VectorArray chain_pos_;

        VectorArray chain_vel_;

        VectorArray chain_acc_;

        IdxArray index_;

        IdxArray new_index_;

        std::conditional_t<Interactions::ext_vel_dep, VectorArray, Empty> aux_vel_;

        std::conditional_t<Interactions::ext_vel_dep, VectorArray, Empty> chain_aux_vel_;
    };

    /*---------------------------------------------------------------------------*\
        Class ARchainSystem Implementation
    \*---------------------------------------------------------------------------*/
    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    template <CONCEPT_PARTICLE_CONTAINER STL>
    ARchainSystem<Particles, Interactions, RegType>::ARchainSystem(Scalar time, const STL &particle_set)
        : ptcl_(time, particle_set),
          accels_(particle_set.size()),
          regu_(ptcl_),
          chain_pos_(particle_set.size()),
          chain_vel_(particle_set.size()),
          chain_acc_(particle_set.size()),
          index_(particle_set.size()),
          new_index_(particle_set.size()) {
        static_assert(is_ranges_v<STL>, "Only STL-like container can be used");
        Chain::calc_chain_index(ptcl_.pos(), index_);
        Chain::calc_chain(ptcl_.pos(), chain_pos(), index());
        Chain::calc_chain(ptcl_.vel(), chain_vel(), index());
        if constexpr (Interactions::ext_vel_dep) {
            aux_vel_ = ptcl_.vel();
            chain_aux_vel_ = chain_vel_;
        }
        regu_ = std::move(
            Regularization<Scalar, RegType>{*this});  // re-construct the regularization with chain coordinates.
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    void ARchainSystem<Particles, Interactions, RegType>::advance_time(Scalar step_size) {
        Scalar phy_time = regu_.eval_pos_phy_time(*this, step_size);
        ptcl_.time() += phy_time;
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    void ARchainSystem<Particles, Interactions, RegType>::advance_pos(Scalar step_size, VectorArray const &velocity) {
        Scalar phy_time = regu_.eval_pos_phy_time(*this, step_size);
        Chain::calc_chain(velocity, chain_acc_, index());  // borrow chain_acc_ as chain velocity buffer.
        chain_advance(ptcl_.pos(), chain_pos(), chain_acc_, phy_time);
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    void ARchainSystem<Particles, Interactions, RegType>::advance_vel(Scalar step_size,
                                                                      VectorArray const &acceleration) {
        Scalar phy_time = regu_.eval_vel_phy_time(*this, step_size);
        Chain::calc_chain(acceleration, chain_acc_, index());
        chain_advance(ptcl_.vel(), chain_vel(), chain_acc_, phy_time);
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    void ARchainSystem<Particles, Interactions, RegType>::evaluate_acc(VectorArray &acceleration) const {
        Interactions::eval_acc(*this, acceleration);
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    void ARchainSystem<Particles, Interactions, RegType>::drift(Scalar step_size) {
        Scalar phy_time = regu_.eval_pos_phy_time(*this, step_size);
        chain_advance(ptcl_.pos(), chain_pos(), chain_vel(), phy_time);
        ptcl_.time() += phy_time;
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
            chain_advance(ptcl_.vel(), chain_vel(), chain_acc_, half_time);
            advance_omega(ptcl_.vel(), accels_.newtonian_acc(), phy_time);
            if constexpr (Interactions::ext_vel_indep) {
                advance_bindE(ptcl_.vel(), accels_.ext_vel_indep_acc(), phy_time);
            }
            chain_advance(ptcl_.vel(), chain_vel(), chain_acc_, half_time);

            /*advance_omega(ptcl_.vel(), accels_.newtonian_acc(), half_time);
            if constexpr (Interactions::ext_vel_indep) {
                advance_bindE(ptcl_.vel(), accels_.ext_vel_indep_acc(), half_time);
            }

            Chain::calc_chain(accels_.tot_vel_indep_acc(), chain_acc_, index());
            chain_advance(ptcl_.vel(), chain_vel(), chain_acc_, phy_time);

            advance_omega(ptcl_.vel(), accels_.newtonian_acc(), half_time);
            if constexpr (Interactions::ext_vel_indep) {
                advance_bindE(ptcl_.vel(), accels_.ext_vel_indep_acc(), half_time);
            }*/
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    void ARchainSystem<Particles, Interactions, RegType>::pre_iter_process() {
        if constexpr (Interactions::ext_vel_dep) {
            aux_vel_ = ptcl_.vel();
            chain_aux_vel_ = chain_vel_;
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    void ARchainSystem<Particles, Interactions, RegType>::post_iter_process() {
        Chain::calc_chain_index(ptcl_.pos(), new_index_);
        if (new_index_ != index_) {
            Chain::update_chain(chain_pos_, index_, new_index_);
            Chain::calc_cartesian(ptcl_.mass(), chain_pos_, ptcl_.pos(), new_index_);
            Chain::update_chain(chain_vel_, index_, new_index_);
            Chain::calc_cartesian(ptcl_.mass(), chain_vel_, ptcl_.vel(), new_index_);
            index_ = new_index_;
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    template <typename STL>
    void ARchainSystem<Particles, Interactions, RegType>::write_to_scalar_array(STL &stl_ranges) {
        stl_ranges.clear();

        stl_ranges.reserve(ptcl_.number() * 6 + 3);
        stl_ranges.emplace_back(ptcl_.time());
        stl_ranges.emplace_back(omega());
        stl_ranges.emplace_back(bindE());
        add_coords_to(stl_ranges, chain_pos_);
        add_coords_to(stl_ranges, chain_vel_);

        /*stl_ranges.reserve(ptcl_.number() * 12 + 3);
        stl_ranges.emplace_back(ptcl_.time());
        stl_ranges.emplace_back(omega());
        stl_ranges.emplace_back(bindE());
        add_coords_to(stl_ranges, chain_pos_);
        add_coords_to(stl_ranges, chain_vel_);
        add_coords_to(stl_ranges, ptcl_.pos());
        add_coords_to(stl_ranges, ptcl_.vel());*/
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    template <typename STL>
    void ARchainSystem<Particles, Interactions, RegType>::read_from_scalar_array(const STL &stl_ranges) {
        auto begin = stl_ranges.begin();
        ptcl_.time() = *begin;
        omega() = *(begin + 1);
        bindE() = *(begin + 2);

        size_t len = ptcl_.number() * 3;

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

        load_to_coords(c_pos_begin, c_pos_end, ptcl_.pos());
        load_to_coords(c_vel_begin, c_vel_end, ptcl_.vel());*/

        Chain::calc_cartesian(ptcl_.mass(), chain_pos_, ptcl_.pos(), index());
        Chain::calc_cartesian(ptcl_.mass(), chain_vel_, ptcl_.vel(), index());
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    std::string ARchainSystem<Particles, Interactions, RegType>::column_names() const {
        return ptcl_.column_names();
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    std::ostream &operator<<(std::ostream &os, ARchainSystem<Particles, Interactions, RegType> const &ps) {
        os << ps.ptcl_;
        return os;
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    std::istream &operator>>(std::istream &is, ARchainSystem<Particles, Interactions, RegType> &ps) {
        is >> ps.ptcl_;
        return is;
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    void ARchainSystem<Particles, Interactions, RegType>::chain_advance(VectorArray &var, VectorArray &chain_var,
                                                                        VectorArray const &chain_increment,
                                                                        Scalar phy_time) {
        calc::array_advance(chain_var, chain_increment, phy_time);
        Chain::calc_cartesian(ptcl_.mass(), chain_var, var, index());
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
    void ARchainSystem<Particles, Interactions, RegType>::advance_omega(VectorArray const &velocity,
                                                                        VectorArray const &d_omega_dr,
                                                                        Scalar phy_time) {
        // if constexpr (regu_type == ReguType::TTL) {
        Scalar d_omega = calc::coord_contract_to_scalar(ptcl_.mass(), velocity, d_omega_dr);
        regu_.omega() += d_omega * phy_time;
        //}
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    void ARchainSystem<Particles, Interactions, RegType>::advance_bindE(VectorArray const &velocity,
                                                                        VectorArray const &d_bindE_dr,
                                                                        Scalar phy_time) {
        if constexpr ((Interactions::ext_vel_indep || Interactions::ext_vel_dep) && regu_type == ReguType::LogH) {
            Scalar d_bindE = -calc::coord_contract_to_scalar(ptcl_.mass(), velocity, d_bindE_dr);
            regu_.bindE() += d_bindE * phy_time;
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    void ARchainSystem<Particles, Interactions, RegType>::kick_pseu_vel(Scalar phy_time) {
        Interactions::eval_extra_vel_dep_acc(*this, accels_.ext_vel_dep_acc());
        calc::array_add(accels_.acc(), accels_.tot_vel_indep_acc(), accels_.ext_vel_dep_acc());
        Chain::calc_chain(accels_.acc(), chain_acc_, index());
        chain_advance(aux_vel_, chain_aux_vel_, chain_acc_, phy_time);

        advance_omega(ptcl_.vel(), accels_.newtonian_acc(), phy_time);

        if constexpr (Interactions::ext_vel_indep) {
            calc::array_add(accels_.acc(), accels_.ext_vel_indep_acc(), accels_.ext_vel_dep_acc());
            advance_bindE(ptcl_.vel(), accels_.acc(), phy_time);
        } else {
            advance_bindE(ptcl_.vel(), accels_.ext_vel_dep_acc(), phy_time);
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    void ARchainSystem<Particles, Interactions, RegType>::kick_real_vel(Scalar phy_time) {
        std::swap(aux_vel_, ptcl_.vel());
        std::swap(chain_aux_vel_, chain_vel());
        Interactions::eval_extra_vel_dep_acc(*this, accels_.ext_vel_dep_acc());
        std::swap(aux_vel_, ptcl_.vel());
        std::swap(chain_aux_vel_, chain_vel());
        calc::array_add(accels_.acc(), accels_.tot_vel_indep_acc(), accels_.ext_vel_dep_acc());

        Chain::calc_chain(accels_.acc(), chain_acc_, index());
        chain_advance(ptcl_.vel(), chain_vel(), chain_acc_, phy_time);
    }
}  // namespace space::particle_system
