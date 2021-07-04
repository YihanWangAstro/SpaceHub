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
 * @file tide-particles.hpp
 *
 * Header file.
 */
#pragma once

#include "../IO.hpp"
#include "../spacehub-concepts.hpp"
#include "../vector/vector3.hpp"
namespace hub::particles {

    /*---------------------------------------------------------------------------*\
    Class TideParticle Declaration
    \*---------------------------------------------------------------------------*/
    /**
     * @brief
     *
     * @tparam Vec3
     */
    template <typename Vec3>
    struct TideParticle {
       public:
        /**
         * @brief 3D vector type.
         */
        using Vector = Vec3;

        /**
         * @brief Floating point like type;
         */
        using Scalar = typename Vector::value_type;

        SPACEHUB_MAKE_CONSTRUCTORS(TideParticle, default, default, default, default, default);

        explicit TideParticle(Scalar mass, Scalar r, Scalar apsidal_motion_const, Scalar lag_time, Vector position,
                              Vector velocity);

        explicit TideParticle(Scalar mass, Scalar r, Scalar apsidal_motion_const, Scalar lag_time, Scalar px = 0,
                              Scalar py = 0, Scalar pz = 0, Scalar vx = 0, Scalar vy = 0, Scalar vz = 0);

        // public members
        /**
         * @brief Position of the particle.
         */
        Vector pos;

        /**
         * @brief Velocity of the particle.
         */
        Vector vel;

        /**
         * @brief Mass of the particle.
         */
        Scalar mass;

        /**
         * @brief
         *
         */
        Scalar radius;

        /**
         * @brief Apsidal motion constant
         */
        Scalar tide_apsidal_const;

        /**
         * @brief Tidal time lag
         */
        Scalar tide_lag_time;
    };

    /*---------------------------------------------------------------------------*\
        Class TideParticles Declaration
    \*---------------------------------------------------------------------------*/
    /**
     * @brief Structure of Array point particle group.
     *
     * @tparam TypeSystem The type system in spaceHub(hub::Types).
     */
    template <typename TypeSystem>
    class TideParticles {
       public:
        // Type members
        SPACEHUB_USING_TYPE_SYSTEM_OF(TypeSystem);

        /**
         * @brief Embedded Particle type.
         */
        using Particle = TideParticle<Vector>;

        // Constructors
        SPACEHUB_MAKE_CONSTRUCTORS(TideParticles, default, default, default, default, default);

        /**
         * @brief Construct a new Point Particles object from std::ranges(Container)
         *
         * @tparam STL std::ranges(Container)
         * @param[in] t Initial time of the the particle group.
         * @param[in] particle_set Input particle set.
         */
        template <CONCEPT_PARTICLE_CONTAINER STL>
        TideParticles(Scalar t, STL const &particle_set);

        // Public methods
        SPACEHUB_STD_ACCESSOR(StateScalar, time, time_);

        SPACEHUB_ARRAY_ACCESSOR(ScalarArray, mass, mass_);

        SPACEHUB_ARRAY_ACCESSOR(ScalarArray, radius, radius_);

        SPACEHUB_ARRAY_ACCESSOR(IdxArray, idn, idn_);

        SPACEHUB_ARRAY_ACCESSOR(StateVectorArray, pos, pos_);

        SPACEHUB_ARRAY_ACCESSOR(StateVectorArray, vel, vel_);

        SPACEHUB_ARRAY_ACCESSOR(ScalarArray, tide_apsidal_const, k_AM_);

        SPACEHUB_ARRAY_ACCESSOR(ScalarArray, tide_lag_time, tau_lag_);

        void resize(size_t new_sz);

        void reserve(size_t new_cap);

        void emplace_back(Particle const &new_particle);

        size_t number() const;

        size_t capacity() const;

        void clear();

        std::string column_names() const;

        std::vector<Particle> to_AoS() const;

        template <typename U>
        friend std::ostream &operator<<(std::ostream &os, TideParticles<U> const &ps);

        template <typename U>
        friend std::istream &operator>>(std::istream &is, TideParticles<U> &ps);

       private:
        // Private members
        StateVectorArray pos_;

        StateVectorArray vel_;

        ScalarArray mass_;

        ScalarArray radius_;

        ScalarArray k_AM_;

        ScalarArray tau_lag_;

        IdxArray idn_;

        StateScalar time_{0.0};

        size_t active_num_{0};
    };
}  // namespace hub::particles

namespace hub::particles {

    /*---------------------------------------------------------------------------*\
        Class TideParticle Implementation
    \*---------------------------------------------------------------------------*/
    template <typename Vec3>
    TideParticle<Vec3>::TideParticle(Scalar m, Scalar r, Scalar apsidal_motion_const, Scalar lag_time, Vec3 position,
                                     Vec3 velocity)
        : pos(position),
          vel(velocity),
          mass(m),
          radius(r),
          tide_apsidal_const(apsidal_motion_const),
          tide_lag_time(lag_time) {}

    template <typename Vec3>
    TideParticle<Vec3>::TideParticle(Scalar m, Scalar r, Scalar apsidal_motion_const, Scalar lag_time, Scalar px,
                                     Scalar py, Scalar pz, Scalar vx, Scalar vy, Scalar vz)
        : pos(px, py, pz),
          vel(vx, vy, vz),
          mass(m),
          radius(r),
          tide_apsidal_const(apsidal_motion_const),
          tide_lag_time(lag_time) {}

    template <typename Vec3>
    std::ostream &operator<<(std::ostream &os, TideParticle<Vec3> const &particle) {
        hub::print_csv(os, particle.mass, particle.radius, particle.pos, particle.vel, particle.tide_apsidal_const,
                       particle.tide_lag_time);
        return os;
    }

    template <typename Vec3>
    std::istream &operator>>(std::istream &is, TideParticle<Vec3> &particle) {
        hub::input(is, particle.mass, particle.radius, particle.pos, particle.vel, particle.tide_apsidal_const,
                   particle.tide_lag_time);
        return is;
    }

    /*---------------------------------------------------------------------------*\
        Class TideParticles Implementation
    \*---------------------------------------------------------------------------*/
    template <typename TypeSystem>
    template <CONCEPT_PARTICLE_CONTAINER STL>
    TideParticles<TypeSystem>::TideParticles(Scalar t, const STL &particle_set) {
        size_t input_num = particle_set.size();
        this->reserve(input_num);
        size_t id = 0;
        for (auto &p : particle_set) {
            pos_.emplace_back(p.pos);
            vel_.emplace_back(p.vel);
            mass_.emplace_back(p.mass);
            radius_.emplace_back(p.radius);
            k_AM_.emplace_back(p.tide_apsidal_const);
            tau_lag_.emplace_back(p.tide_lag_time);
            idn_.emplace_back(id++);
        }
        time_ = t;
        active_num_ = input_num;
    }

    template <typename TypeSystem>
    void TideParticles<TypeSystem>::resize(size_t new_sz) {
        hub::resize_all(new_sz, pos_, vel_, mass_, radius_, k_AM_, tau_lag_, idn_);
        active_num_ = new_sz;
    }

    template <typename TypeSystem>
    void TideParticles<TypeSystem>::reserve(size_t new_cap) {
        hub::reserve_all(new_cap, pos_, vel_, mass_, radius_, k_AM_, tau_lag_, idn_);
    }

    template <typename TypeSystem>
    void TideParticles<TypeSystem>::clear() {
        hub::clear_all(pos_, vel_, mass_, radius_, k_AM_, tau_lag_, idn_);
        active_num_ = 0;
    }

    template <typename TypeSystem>
    void TideParticles<TypeSystem>::emplace_back(typename TideParticles<TypeSystem>::Particle const &new_particle) {
        pos_.emplace_back(new_particle.pos);
        vel_.emplace_back(new_particle.vel);
        mass_.emplace_back(new_particle.mass);
        radius_.emplace_back(new_particle.radius);
        k_AM_.emplace_back(new_particle.tide_apsidal_const);
        tau_lag_.emplace_back(new_particle.tide_lag_time);
        idn_.emplace_back(this->number());
        active_num_++;
    }

    template <typename TypeSystem>
    size_t TideParticles<TypeSystem>::number() const {
        return active_num_;
    }

    template <typename TypeSystem>
    size_t TideParticles<TypeSystem>::capacity() const {
        return idn_.capacity();
    }

    template <typename TypeSystem>
    std::string TideParticles<TypeSystem>::column_names() const {
        return "time,id,mass,radius,k_AM,tau_lag,px,py,pz,vx,vy,vz";
    }

    template <typename TypeSystem>
    auto TideParticles<TypeSystem>::to_AoS() const -> std::vector<Particle> {
        std::vector<Particle> ptc;
        size_t ptc_num = this->number();
        ptc.reserve(ptc_num);
        for (size_t i = 0; i < ptc_num; ++i) {
            ptc.emplace_back(mass_[i], radius_[i], k_AM_[i], tau_lag_[i], pos_[i], vel_[i]);
        }
        return ptc;
    }

    template <typename TypeSystem>
    std::ostream &operator<<(std::ostream &os, TideParticles<TypeSystem> const &ps) {
        size_t num = ps.number();
        for (size_t i = 0; i < num; ++i) {
            hub::print(os, ps.time(), ',', ps.idn(i), ',', ps.mass(i), ',', ps.radius(i), ',', ps.tide_apsidal_const(i),
                       ',', ps.tide_lag_time(i), ',', ps.pos(i), ',', ps.vel(i), '\n');
        }
        return os;
    }
}  // namespace hub::particles