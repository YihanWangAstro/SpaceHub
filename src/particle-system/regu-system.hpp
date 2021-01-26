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
 * @file regu-system.hpp
 *
 * Header file.
 */
#pragma once

#include <type_traits>

#include "../core-computation.hpp"
#include "../interaction/interaction.hpp"
#include "../spacehub-concepts.hpp"
namespace space::particle_system {
    /**
     *
     */
    enum class ReguType { LogH, TTL, None };

    /*---------------------------------------------------------------------------*\
        Class Regularization Declaration
    \*---------------------------------------------------------------------------*/
    /**
     *
     * @tparam
     * @tparam Type
     */
    template <typename TypeSystem, ReguType Type = ReguType::LogH>
    class Regularization {
       public:
        SPACEHUB_USING_TYPE_SYSTEM_OF(TypeSystem);

        // Constructors
        template <typename Particles>
        explicit Regularization(Particles const &particles);

        // Public methods
        SPACEHUB_STD_ACCESSOR(Scalar, omega, omega_);

        SPACEHUB_STD_ACCESSOR(Scalar, bindE, bindE_);

        template <CONCEPT_PARTICLES_DATA Particles>
        Scalar eval_pos_phy_time(Particles const &particles, Scalar step_size) const;

        template <CONCEPT_PARTICLES_DATA Particles>
        Scalar eval_vel_phy_time(Particles const &particles, Scalar step_size) const;

        template <typename Particles>
        inline Scalar regu_function(Particles const &particles) const;

       private:
        template <typename Particles>
        inline Scalar capital_omega(Particles const &particles) const;
        // Private members
        Scalar omega_{1};

        Scalar bindE_{1};

        Scalar scale_{1};
    };

    /*---------------------------------------------------------------------------*\
        Class RegularizedSystem Declaration
    \*---------------------------------------------------------------------------*/

    /**
     * @brief Regularized particle System.
     *
     * Regularied particle system.
     * See details [Explicit Symplectic Algorithms For Time‐Transformed
     * Hamiltonians](https://link.springer.com/article/10.1023%2FA%3A1008368322547) , [A Class of Symplectic Integrators
     * with Adaptive Time Step for Separable Hamiltonian Systems](http://iopscience.iop.org/article/10.1086/301102/meta)
     * and [A Time-Transformed Leapfrog Scheme](https://link.springer.com/article/10.1023%2FA%3A1021149313347).
     * @tparam Particles
     * @tparam Interactions
     */
    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    class RegularizedSystem : public Particles {
       public:
        // Type members
        SPACEHUB_USING_TYPE_SYSTEM_OF(Particles);

        using Particle = typename Particles::Particle;

        using Interaction = Interactions;

        // Constructors
        SPACEHUB_MAKE_CONSTRUCTORS(RegularizedSystem, delete, default, default, default, default);

        template <CONCEPT_PARTICLE_CONTAINER STL>
        RegularizedSystem(Scalar time, STL const &particle_set);

        // Static members
        static constexpr ReguType regu_type{RegType};

        // Public methods
        SPACEHUB_STD_ACCESSOR(auto, omega, regu_.omega());

        SPACEHUB_STD_ACCESSOR(auto, bindE, regu_.bindE());

        Scalar regu_function() const { return regu_.regu_function(*this); };

        void advance_time(Scalar step_size);

        void advance_pos(Scalar step_size, VectorArray const &velocity);

        void advance_vel(Scalar step_size, VectorArray const &acceleration);

        void evaluate_acc(VectorArray &acceleration) const;

        void drift(Scalar step_size);

        void kick(Scalar step_size);

        void pre_iter_process();

        void post_iter_process(){};

        template <typename STL>
        void write_to_scalar_array(STL &stl_ranges);

        template <typename STL>
        void read_from_scalar_array(STL const &stl_ranges);

        /**
         * @brief
         *
         * @tparam STL
         * @param stl_ranges
         */
        template <typename STL>
        void evaluate_general_derivative(STL &stl_ranges) const;

        // Friend functions
        template <CONCEPT_PARTICLES P, CONCEPT_INTERACTION F, ReguType R>
        friend std::ostream &operator<<(std::ostream &os, RegularizedSystem<P, F, R> const &ps);

        template <CONCEPT_PARTICLES P, CONCEPT_INTERACTION F, ReguType R>
        friend std::istream &operator>>(std::istream &is, RegularizedSystem<P, F, R> &ps);

       private:
        // Private methods
        void eval_vel_indep_acc();

        void advance_omega(VectorArray const &velocity, VectorArray const &d_omega_dr, Scalar phy_time);

        void advance_bindE(VectorArray const &velocity, VectorArray const &d_bindE_dr, Scalar phy_time);

        void kick_pseu_vel(Scalar phy_time);

        void kick_real_vel(Scalar phy_time);

        // Private members
        interactions::InteractionData<Interactions, VectorArray> accels_;
        Regularization<TypeSet, RegType> regu_;
        std::conditional_t<Interactions::ext_vel_dep, VectorArray, Empty> aux_vel_;
    };

    /*---------------------------------------------------------------------------*\
        Class RegularizedSystem Implementation
    \*---------------------------------------------------------------------------*/
    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    template <CONCEPT_PARTICLE_CONTAINER STL>
    RegularizedSystem<Particles, Interactions, RegType>::RegularizedSystem(Scalar time, const STL &particle_set)
        : Particles(time, particle_set), accels_(particle_set.size()), regu_(*this) {
        if constexpr (Interactions::ext_vel_dep) {
            aux_vel_ = this->vel();
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    std::istream &operator>>(std::istream &is, RegularizedSystem<Particles, Interactions, RegType> &ps) {
        is >> static_cast<Particles>(ps);
        return is;
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    std::ostream &operator<<(std::ostream &os, const RegularizedSystem<Particles, Interactions, RegType> &ps) {
        os << static_cast<Particles>(ps);
        return os;
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    void RegularizedSystem<Particles, Interactions, RegType>::advance_time(Scalar step_size) {
        Scalar phy_time = regu_.eval_pos_phy_time(*this, step_size);
        this->time() += phy_time;
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    void RegularizedSystem<Particles, Interactions, RegType>::advance_pos(Scalar step_size,
                                                                          VectorArray const &velocity) {
        Scalar phy_time = regu_.eval_pos_phy_time(*this, step_size);
        calc::array_advance(this->pos(), velocity, phy_time);
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    void RegularizedSystem<Particles, Interactions, RegType>::advance_vel(Scalar step_size,
                                                                          VectorArray const &acceleration) {
        Scalar phy_time = regu_.eval_vel_phy_time(*this, step_size);
        Scalar half_time = 0.5 * phy_time;

        if constexpr (regu_type == ReguType::TTL) {
            calc::array_advance(this->vel(), acceleration, half_time);
            Interactions::eval_newtonian_acc(*this, accels_.newtonian_acc());
            advance_omega(this->vel(), accels_.newtonian_acc(), phy_time);
            calc::array_advance(this->vel(), acceleration, half_time);
        } else if constexpr (regu_type == ReguType::LogH) {
            if constexpr (Interactions::ext_vel_indep || Interactions::ext_vel_dep) {
                calc::array_advance(this->vel(), acceleration, half_time);
                Interactions::eval_extra_acc(*this, accels_.acc());
                advance_bindE(this->vel(), accels_.acc(), phy_time);
                calc::array_advance(this->vel(), acceleration, half_time);
            } else {
                calc::array_advance(this->vel(), acceleration, phy_time);
            }
        } else {
            calc::array_advance(this->vel(), acceleration, phy_time);
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    void RegularizedSystem<Particles, Interactions, RegType>::evaluate_acc(VectorArray &acceleration) const {
        Interactions::eval_acc(*this, acceleration);
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    void RegularizedSystem<Particles, Interactions, RegType>::drift(Scalar step_size) {
        Scalar phy_time = regu_.eval_pos_phy_time(*this, step_size);
        calc::array_advance(this->pos(), this->vel(), phy_time);
        this->time() += phy_time;
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    void RegularizedSystem<Particles, Interactions, RegType>::kick(Scalar step_size) {
        Scalar phy_time = regu_.eval_vel_phy_time(*this, step_size);
        Scalar half_time = 0.5 * phy_time;

        eval_vel_indep_acc();

        if constexpr (Interactions::ext_vel_dep) {
            kick_real_vel(half_time);
            kick_pseu_vel(phy_time);
            kick_real_vel(half_time);
        } else {
            advance_omega(this->vel(), accels_.newtonian_acc(), half_time);
            if constexpr (Interactions::ext_vel_indep) {
                advance_bindE(this->vel(), accels_.ext_vel_indep_acc(), half_time);
            }
            calc::array_advance(this->vel(), accels_.tot_vel_indep_acc(), phy_time);
            if constexpr (Interactions::ext_vel_indep) {
                advance_bindE(this->vel(), accels_.ext_vel_indep_acc(), half_time);
            }
            advance_omega(this->vel(), accels_.newtonian_acc(), half_time);
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    void RegularizedSystem<Particles, Interactions, RegType>::pre_iter_process() {
        if constexpr (Interactions::ext_vel_dep) {
            aux_vel_ = this->vel();
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    template <typename STL>
    void RegularizedSystem<Particles, Interactions, RegType>::write_to_scalar_array(STL &stl_ranges) {
        stl_ranges.clear();
        stl_ranges.reserve(this->number() * 3 * (2 + static_cast<size_t>(Interactions::ext_vel_dep)) + 3);
        stl_ranges.emplace_back(this->time());
        stl_ranges.emplace_back(omega());
        stl_ranges.emplace_back(bindE());

        add_coords_to(stl_ranges, this->pos());
        add_coords_to(stl_ranges, this->vel());
        if constexpr (Interactions::ext_vel_dep) {
            add_coords_to(stl_ranges, aux_vel_);
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    template <typename STL>
    void RegularizedSystem<Particles, Interactions, RegType>::evaluate_general_derivative(STL &stl_ranges) const {
        stl_ranges.clear();
        stl_ranges.reserve(this->number() * 3 * (2 + static_cast<size_t>(Interactions::ext_vel_dep)) + 3);

        Scalar pos_regu = regu_.eval_pos_phy_time(1);
        Scalar vel_regu = regu_.eval_vel_phy_time(1);

        stl_ranges.emplace_back(pos_regu);

        Interactions::eval_newtonian_acc(*this, accels_.newtonian_acc());

        if constexpr (regu_type == ReguType::TTL) {
            Scalar d_omega_dh =
                calc::coord_contract_to_scalar(this->mass(), this->vel(), accels_.newtonian_acc()) * vel_regu;
            stl_ranges.emplace_back(d_omega_dh);
        } else {
            stl_ranges.emplace_back(0);
        }

        if constexpr ((Interactions::ext_vel_indep || Interactions::ext_vel_dep) && regu_type == ReguType::LogH) {
            Interactions::eval_extra_acc(*this, accels_.acc());
            Scalar d_bindE_dh = -calc::coord_contract_to_scalar(this->mass(), this->vel(), accels_.acc()) * vel_regu;
            stl_ranges.emplace_back(d_bindE_dh);
            calc::array_add(accels_.acc(), accels_.acc(), accels_.newtonian_acc());
        } else {
            stl_ranges.emplace_back(0);
        }

        add_scaled_coords_to(stl_ranges, this->vel(), pos_regu);
        add_scaled_coords_to(stl_ranges, accels_.acc(), vel_regu);
        if constexpr (Interactions::ext_vel_dep) {
            add_scaled_coords_to(stl_ranges, accels_.acc(), vel_regu);
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    template <typename STL>
    void RegularizedSystem<Particles, Interactions, RegType>::read_from_scalar_array(const STL &stl_ranges) {
        auto begin = stl_ranges.begin();
        this->time() = *begin;
        omega() = *(++begin);
        bindE() = *(++begin);
        size_t len = this->number() * 3;
        auto pos_begin = ++begin;
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

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    void RegularizedSystem<Particles, Interactions, RegType>::eval_vel_indep_acc() {
        Interactions::eval_newtonian_acc(*this, accels_.newtonian_acc());

        if constexpr (Interactions::ext_vel_indep) {
            Interactions::eval_extra_vel_indep_acc(*this, accels_.ext_vel_indep_acc());
            calc::array_add(accels_.tot_vel_indep_acc(), accels_.ext_vel_indep_acc(), accels_.newtonian_acc());
        } else {
            accels_.tot_vel_indep_acc() = accels_.newtonian_acc();
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    void RegularizedSystem<Particles, Interactions, RegType>::advance_omega(VectorArray const &velocity,
                                                                            VectorArray const &d_omega_dr,
                                                                            Scalar phy_time) {
        if constexpr (regu_type == ReguType::TTL) {
            Scalar d_omega = calc::coord_contract_to_scalar(this->mass(), velocity, d_omega_dr);
            regu_.omega() += d_omega * phy_time;
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    void RegularizedSystem<Particles, Interactions, RegType>::advance_bindE(VectorArray const &velocity,
                                                                            VectorArray const &d_bindE_dr,
                                                                            Scalar phy_time) {
        if constexpr ((Interactions::ext_vel_indep || Interactions::ext_vel_dep) && regu_type == ReguType::LogH) {
            Scalar d_bindE = -calc::coord_contract_to_scalar(this->mass(), velocity, d_bindE_dr);
            regu_.bindE() += d_bindE * phy_time;
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    void RegularizedSystem<Particles, Interactions, RegType>::kick_pseu_vel(Scalar phy_time) {
        Interactions::eval_extra_vel_dep_acc(*this, accels_.ext_vel_dep_acc());
        calc::array_add(accels_.acc(), accels_.tot_vel_indep_acc(), accels_.ext_vel_dep_acc());
        calc::array_advance(aux_vel_, accels_.acc(), phy_time);

        advance_omega(this->vel(), accels_.newtonian_acc(), phy_time);

        if constexpr (Interactions::ext_vel_indep) {
            calc::array_add(accels_.acc(), accels_.ext_vel_indep_acc(), accels_.ext_vel_dep_acc());
            advance_bindE(this->vel(), accels_.acc(), phy_time);
        } else {
            advance_bindE(this->vel(), accels_.ext_vel_dep_acc(), phy_time);
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    void RegularizedSystem<Particles, Interactions, RegType>::kick_real_vel(Scalar phy_time) {
        std::swap(aux_vel_, this->vel());
        Interactions::eval_extra_vel_dep_acc(*this, accels_.ext_vel_dep_acc());
        std::swap(aux_vel_, this->vel());

        calc::array_add(accels_.acc(), accels_.tot_vel_indep_acc(), accels_.ext_vel_dep_acc());
        calc::array_advance(this->vel(), accels_.acc(), phy_time);
    }

    /*---------------------------------------------------------------------------*\
        Class Regularization Implementation
    \*---------------------------------------------------------------------------*/
    template <typename TypeSystem, ReguType Type>
    template <typename Particles>
    Regularization<TypeSystem, Type>::Regularization(Particles const &particles) {
        omega_ = capital_omega(particles);
        bindE_ = -calc::calc_total_energy(particles);
    }

    template <typename TypeSystem, ReguType Type>
    template <CONCEPT_PARTICLES_DATA Particles>
    auto Regularization<TypeSystem, Type>::eval_pos_phy_time(Particles const &particles, Scalar step_size) const
        -> Scalar {
        if constexpr (Type == ReguType::LogH) {
            return step_size / (bindE_ + calc::calc_kinetic_energy(particles));
        } else if constexpr (Type == ReguType::TTL) {
            return step_size / omega_;
        } else if constexpr (Type == ReguType::None) {
            return step_size;
        } else {
            spacehub_abort("Undefined regularization type!");
        }
    }

    template <typename TypeSystem, ReguType Type>
    template <CONCEPT_PARTICLES_DATA Particles>
    auto Regularization<TypeSystem, Type>::eval_vel_phy_time(Particles const &particles, Scalar step_size) const
        -> Scalar {
        if constexpr (Type == ReguType::LogH) {
            return step_size / -calc::calc_potential_energy(particles);
        } else if constexpr (Type == ReguType::TTL) {
            return step_size / capital_omega(particles);
        } else if constexpr (Type == ReguType::None) {
            return step_size;
        } else {
            spacehub_abort("Undefined regularization type!");
        }
    }

    template <typename TypeSystem, ReguType Type>
    template <typename Particles>
    auto Regularization<TypeSystem, Type>::capital_omega(Particles const &particles) const -> Scalar {
        return -calc::calc_potential_energy(particles);
    }

    template <typename TypeSystem, ReguType Type>
    template <typename Particles>
    auto Regularization<TypeSystem, Type>::regu_function(Particles const &particles) const -> Scalar {
        return capital_omega(particles);
    }
}  // namespace space::particle_system
