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
#include "../type-class.hpp"
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
        SPACEHUB_STD_ACCESSOR(StateScalar, omega, omega_);

        SPACEHUB_STD_ACCESSOR(StateScalar, bindE, bindE_);

        template <CONCEPT_PARTICLES_DATA Particles>
        Scalar eval_pos_phy_time(Particles const &particles, Scalar step_size) const;

        template <CONCEPT_PARTICLES_DATA Particles>
        Scalar eval_vel_phy_time(Particles const &particles, Scalar step_size) const;

        template <typename Particles>
        inline StateScalar regu_function(Particles const &particles) const;

       private:
        template <typename Particles>
        inline StateScalar capital_omega(Particles const &particles) const;

        // Private members
        StateScalar omega_{1};

        StateScalar bindE_{1};

        Scalar scale_{1};
    };

    /*---------------------------------------------------------------------------*\
        Class RegularizedSystem Declaration
    \*---------------------------------------------------------------------------*/

    /**
     * @brief Regularized particle System.
     *
     * Regularized particle system.
     * See details [Explicit Symplectic Algorithms For Time‐Transformed
     * Hamiltonian](https://link.springer.com/article/10.1023%2FA%3A1008368322547) , [A Class of Symplectic Integrators
     * with Adaptive Time Step for Separable Hamiltonian Systems](http://iopscience.iop.org/article/10.1086/301102/meta)
     * and [A Time-Transformed Leapfrog Scheme](https://link.springer.com/article/10.1023%2FA%3A1021149313347).
     * @tparam Particles
     * @tparam Interactions
     */
    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType = ReguType::LogH>
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
        SPACEHUB_STD_ACCESSOR(StateScalar, omega, regu_.omega());

        SPACEHUB_STD_ACCESSOR(StateScalar, bindE, regu_.bindE());

        SPACEHUB_STD_ACCESSOR(StateScalarArray, increment, increment_);

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

        void post_iter_process(){};

        template <typename ScalarIterable>
        void write_to_scalar_array(ScalarIterable &stl_ranges);

        template <typename ScalarIterable>
        void read_from_scalar_array(ScalarIterable const &stl_ranges);

        template <typename ScalarIterable>
        void evaluate_general_derivative(ScalarIterable &stl_ranges);

        inline void collect_increment(bool sync) { sync_increment_ = sync; };

        void clear_increment() { calc::array_set_zero(increment_); };

        [[nodiscard]] size_t variable_number() const;

        // Friend functions
        template <CONCEPT_PARTICLES P, CONCEPT_INTERACTION F, ReguType R>
        friend std::ostream &operator<<(std::ostream &os, RegularizedSystem<P, F, R> const &ps);

        template <CONCEPT_PARTICLES P, CONCEPT_INTERACTION F, ReguType R>
        friend std::istream &operator>>(std::istream &is, RegularizedSystem<P, F, R> &ps);

       private:
        // Private methods
        void eval_vel_indep_acc();

        void advance_omega(StateVectorArray const &velocity, VectorArray const &d_omega_dr, Scalar phy_time);

        void advance_bindE(StateVectorArray const &velocity, VectorArray const &d_bindE_dr, Scalar phy_time);

        void kick_pseu_vel(Scalar phy_time);

        void kick_real_vel(Scalar phy_time);

        template <typename Array>
        void sync_pos_increment(Array const &inc, Scalar step_size);

        template <typename Array>
        void sync_vel_increment(Array const &inc, Scalar step_size);

        template <typename Array>
        void sync_auxi_vel_increment(Array const &inc, Scalar step_size);

        void sync_time_increment(Scalar phy_time);

        void sync_omega_increment(Scalar domega);

        void sync_bindE_increment(Scalar dbindE);

        inline constexpr size_t time_offset() { return 0; };
        inline constexpr size_t omega_offset() { return 1; };
        inline constexpr size_t bindE_offset() { return 2; };
        inline constexpr size_t pos_offset() { return 3; };
        inline constexpr size_t vel_offset() { return this->number() * 3 + 3; };
        inline constexpr size_t auxi_vel_offset() { return this->number() * 6 + 3; };

        // Private members
        force::InteractionData<Interactions, VectorArray> accels_;
        StateScalarArray increment_;
        Regularization<TypeSet, RegType> regu_;
        std::conditional_t<Interactions::ext_vel_dep, StateVectorArray, Empty> aux_vel_;
        bool sync_increment_{false};
    };

    /*---------------------------------------------------------------------------*\
        Class RegularizedSystem Implementation
    \*---------------------------------------------------------------------------*/
    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    template <CONCEPT_PARTICLE_CONTAINER STL>
    RegularizedSystem<Particles, Interactions, RegType>::RegularizedSystem(Scalar time, const STL &particle_set)
        : Particles(time, particle_set),
          accels_(particle_set.size()),
          increment_(this->variable_number()),
          regu_(*this) {
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
        sync_time_increment(phy_time);
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    template <typename GenVectorArray>
    void RegularizedSystem<Particles, Interactions, RegType>::advance_pos(Scalar step_size,
                                                                          GenVectorArray const &velocity) {
        Scalar phy_time = regu_.eval_pos_phy_time(*this, step_size);
        calc::array_advance(this->pos(), velocity, phy_time);
        sync_pos_increment(velocity, phy_time);
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    template <typename GenVectorArray>
    void RegularizedSystem<Particles, Interactions, RegType>::advance_vel(Scalar step_size,
                                                                          GenVectorArray const &acceleration) {
        Scalar phy_time = regu_.eval_vel_phy_time(*this, step_size);
        Scalar half_time = 0.5 * phy_time;

        if constexpr (regu_type == ReguType::TTL) {
            Interactions::eval_newtonian_acc(*this, accels_.newtonian_acc());
            advance_omega(this->vel(), accels_.newtonian_acc(), half_time);
            calc::array_advance(this->vel(), acceleration, phy_time);
            sync_vel_increment(acceleration, phy_time);
            advance_omega(this->vel(), accels_.newtonian_acc(), half_time);
        } else if constexpr (regu_type == ReguType::LogH) {
            if constexpr (Interactions::ext_vel_indep || Interactions::ext_vel_dep) {
                calc::array_advance(this->vel(), acceleration, half_time);
                sync_vel_increment(acceleration, half_time);
                Interactions::eval_extra_acc(*this, accels_.acc());
                advance_bindE(this->vel(), accels_.acc(), phy_time);
                calc::array_advance(this->vel(), acceleration, half_time);
                sync_vel_increment(acceleration, half_time);
            } else {
                calc::array_advance(this->vel(), acceleration, phy_time);
                sync_vel_increment(acceleration, phy_time);
            }
        } else {
            calc::array_advance(this->vel(), acceleration, phy_time);
            sync_vel_increment(acceleration, phy_time);
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    template <typename GenVectorArray>
    void RegularizedSystem<Particles, Interactions, RegType>::evaluate_acc(GenVectorArray &acceleration) const {
        Interactions::eval_acc(*this, acceleration);
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    void RegularizedSystem<Particles, Interactions, RegType>::drift(Scalar step_size) {
        Scalar phy_time = regu_.eval_pos_phy_time(*this, step_size);
        calc::array_advance(this->pos(), this->vel(), phy_time);
        this->time() += phy_time;
        sync_time_increment(phy_time);
        sync_pos_increment(this->vel(), phy_time);
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
            sync_vel_increment(accels_.tot_vel_indep_acc(), phy_time);
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
    template <typename ScalarIterable>
    void RegularizedSystem<Particles, Interactions, RegType>::write_to_scalar_array(ScalarIterable &stl_ranges) {
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
    template <typename ScalarIterable>
    void RegularizedSystem<Particles, Interactions, RegType>::evaluate_general_derivative(ScalarIterable &stl_ranges) {
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

        add_scaled_coords_to(stl_ranges, this->vel(), pos_regu);
        if constexpr (Interactions::ext_vel_indep || Interactions::ext_vel_dep) {
            add_scaled_coords_to(stl_ranges, accels_.acc(), vel_regu);
        } else {
            add_scaled_coords_to(stl_ranges, accels_.newtonian_acc(), vel_regu);
        }
        if constexpr (Interactions::ext_vel_dep) {
            add_scaled_coords_to(stl_ranges, accels_.acc(), vel_regu);
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    template <typename ScalarIterable>
    void RegularizedSystem<Particles, Interactions, RegType>::read_from_scalar_array(const ScalarIterable &stl_ranges) {
        if (stl_ranges.size() == this->variable_number()) {
            auto begin = stl_ranges.begin();
            this->time() = *begin;
            omega() = *(begin + omega_offset());
            bindE() = *(begin + bindE_offset());
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

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    template <typename Array>
    void RegularizedSystem<Particles, Interactions, RegType>::sync_pos_increment(Array const &inc, Scalar step_size) {
        if (sync_increment_) {
            size_t offset = pos_offset();
            for (size_t i = 0; i < inc.size(); ++i) {
                increment_[3 * i + offset] += inc[i].x * step_size;
                increment_[3 * i + offset + 1] += inc[i].y * step_size;
                increment_[3 * i + offset + 2] += inc[i].z * step_size;
            }
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    template <typename Array>
    void RegularizedSystem<Particles, Interactions, RegType>::sync_vel_increment(Array const &inc, Scalar step_size) {
        if (sync_increment_) {
            size_t offset = vel_offset();
            for (size_t i = 0; i < inc.size(); ++i) {
                increment_[3 * i + offset] += inc[i].x * step_size;
                increment_[3 * i + offset + 1] += inc[i].y * step_size;
                increment_[3 * i + offset + 2] += inc[i].z * step_size;
            }
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    template <typename Array>
    void RegularizedSystem<Particles, Interactions, RegType>::sync_auxi_vel_increment(Array const &inc,
                                                                                      Scalar step_size) {
        if (sync_increment_) {
            size_t offset = auxi_vel_offset();
            for (size_t i = 0; i < inc.size(); ++i) {
                increment_[3 * i + offset] += inc[i].x * step_size;
                increment_[3 * i + offset + 1] += inc[i].y * step_size;
                increment_[3 * i + offset + 2] += inc[i].z * step_size;
            }
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    void RegularizedSystem<Particles, Interactions, RegType>::sync_time_increment(Scalar phy_time) {
        if (sync_increment_) {
            increment_[time_offset()] += phy_time;
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    void RegularizedSystem<Particles, Interactions, RegType>::sync_omega_increment(Scalar domega) {
        if (sync_increment_) {
            increment_[omega_offset()] += domega;
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    void RegularizedSystem<Particles, Interactions, RegType>::sync_bindE_increment(Scalar dbindE) {
        if (sync_increment_) {
            increment_[bindE_offset()] += dbindE;
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    size_t RegularizedSystem<Particles, Interactions, RegType>::variable_number() const {
        return this->number() * 3 * (2 + static_cast<size_t>(Interactions::ext_vel_dep)) + 3;
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
    void RegularizedSystem<Particles, Interactions, RegType>::advance_omega(StateVectorArray const &velocity,
                                                                            VectorArray const &d_omega_dr,
                                                                            Scalar phy_time) {
        if constexpr (regu_type == ReguType::TTL) {
            Scalar d_omega = calc::coord_contract_to_scalar(this->mass(), velocity, d_omega_dr) * phy_time;
            regu_.omega() += d_omega;
            sync_omega_increment(d_omega);
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    void RegularizedSystem<Particles, Interactions, RegType>::advance_bindE(StateVectorArray const &velocity,
                                                                            VectorArray const &d_bindE_dr,
                                                                            Scalar phy_time) {
        if constexpr ((Interactions::ext_vel_indep || Interactions::ext_vel_dep) && regu_type == ReguType::LogH) {
            Scalar d_bindE = -calc::coord_contract_to_scalar(this->mass(), velocity, d_bindE_dr) * phy_time;
            regu_.bindE() += d_bindE;
            sync_bindE_increment(d_bindE);
        }
    }

    template <CONCEPT_PARTICLES Particles, CONCEPT_INTERACTION Interactions, ReguType RegType>
    void RegularizedSystem<Particles, Interactions, RegType>::kick_pseu_vel(Scalar phy_time) {
        Interactions::eval_extra_vel_dep_acc(*this, accels_.ext_vel_dep_acc());
        calc::array_add(accels_.acc(), accels_.tot_vel_indep_acc(), accels_.ext_vel_dep_acc());
        calc::array_advance(aux_vel_, accels_.acc(), phy_time);
        sync_auxi_vel_increment(accels_.acc(), phy_time);
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
        sync_vel_increment(accels_.acc(), phy_time);
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
    auto Regularization<TypeSystem, Type>::capital_omega(Particles const &particles) const -> StateScalar {
        return -calc::calc_potential_energy(particles);
    }

    template <typename TypeSystem, ReguType Type>
    template <typename Particles>
    auto Regularization<TypeSystem, Type>::regu_function(Particles const &particles) const -> StateScalar {
        return capital_omega(particles);
    }

}  // namespace space::particle_system
