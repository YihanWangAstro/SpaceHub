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
 * @file point-particles.hpp
 *
 * Header file.
 */
#pragma once

#include "../IO.hpp"
#include "../spacehub-concepts.hpp"
#include "../vector/vector3.hpp"
namespace space::particle_set {

    /*---------------------------------------------------------------------------*\
    Class PointParticle Declaration
    \*---------------------------------------------------------------------------*/
    /**
     * @brief Point particle
     * @tparam Real Floating point like type.
     */
    template <typename Vec3>
    struct PointParticle {
       public:
        /**
         * @brief 3D vector type.
         */
        using Vector = Vec3;

        /**
         * @brief Floating point like type;
         */
        using Scalar = typename Vector::value_type;

        SPACEHUB_MAKE_CONSTRUCTORS(PointParticle, default, default, default, default, default);

        /**
         * @brief Construct a new Point Particle object
         *
         * @param[in] mass The mass of the particle
         * @param[in] position The 3D vector position of the particle
         * @param[in] velocity The 3D vector velocity of the particle
         */
        explicit PointParticle(Scalar mass, Vector position, Vector velocity);

        /**
         * @brief Construct a new Point Particle object
         *
         * @param[in] mass The mass fof the particle
         * @param[in] px The x-component of the position vector
         * @param[in] py The y-component of the position vector
         * @param[in] pz The z-component of the position vector
         * @param[in] vx The x-component of the velocity vector
         * @param[in] vy The y-component of the velocity vector
         * @param[in] vz The z-component of the velocity vector
         */
        explicit PointParticle(Scalar mass, Scalar px = 0, Scalar py = 0, Scalar pz = 0, Scalar vx = 0, Scalar vy = 0,
                               Scalar vz = 0);

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
    };

    /*---------------------------------------------------------------------------*\
        Class PointParticles Declaration
    \*---------------------------------------------------------------------------*/
    /**
     * @brief Structure of Array point particle group.
     *
     * @tparam TypeSystem The type system in spaceHub(space::Types).
     */
    template <typename TypeSystem>
    class PointParticles {
       public:
        // Type members
        SPACEHUB_USING_TYPE_SYSTEM_OF(TypeSystem);

        /**
         * @brief Embedded Particle type.
         */
        using Particle = PointParticle<Vector>;

        // Constructors
        SPACEHUB_MAKE_CONSTRUCTORS(PointParticles, default, default, default, default, default);

        /**
         * @brief Construct a new Point Particles object from std::ranges(Container)
         *
         * @tparam STL std::ranges(Container)
         * @param[in] t Initial time of the the particle group.
         * @param[in] particle_set Input particle set.
         */
        template <CONCEPT_PARTICLE_CONTAINER STL>
        PointParticles(Scalar t, STL const &particle_set);

        // Public methods
        SPACEHUB_STD_ACCESSOR(AdScalar, time, time_);

        SPACEHUB_ARRAY_ACCESSOR(ScalarArray, mass, mass_);

        SPACEHUB_ARRAY_ACCESSOR(IdxArray, idn, idn_);

        SPACEHUB_ARRAY_ACCESSOR(AdVectorArray, pos, pos_);

        SPACEHUB_ARRAY_ACCESSOR(AdVectorArray, vel, vel_);

        void resize(size_t new_sz);

        void reserve(size_t new_cap);

        void emplace_back(Particle const &new_particle);

        [[nodiscard]] size_t number() const;

        [[nodiscard]] size_t capacity() const;

        void clear();

        std::string column_names() const;

        template <typename U>
        friend std::ostream &operator<<(std::ostream &os, PointParticles<U> const &ps);

        template <typename U>
        friend std::istream &operator>>(std::istream &is, PointParticles<U> &ps);

       private:
        // Private members
        AdVectorArray pos_;

        AdVectorArray vel_;

        ScalarArray mass_;

        IdxArray idn_;

        AdScalar time_{0.0};

        size_t active_num_{0};
    };
}  // namespace space::particle_set

namespace space::particle_set {

    /*---------------------------------------------------------------------------*\
        Class PointParticle Implementation
    \*---------------------------------------------------------------------------*/
    template <typename Vec3>
    PointParticle<Vec3>::PointParticle(Scalar m, Vec3 position, Vec3 velocity)
        : pos(position), vel(velocity), mass(m) {}

    template <typename Vec3>
    PointParticle<Vec3>::PointParticle(Scalar m, Scalar px, Scalar py, Scalar pz, Scalar vx, Scalar vy, Scalar vz)
        : pos(px, py, pz), vel(vx, vy, vz), mass(m) {}

    template <typename Vec3>
    std::ostream &operator<<(std::ostream &os, PointParticle<Vec3> const &particle) {
        space::print_csv(os, particle.mass, particle.pos, particle.vel);
        return os;
    }

    template <typename Vec3>
    std::istream &operator>>(std::istream &is, PointParticle<Vec3> &particle) {
        space::input(is, particle.mass, particle.pos, particle.vel);
        return is;
    }

    /*---------------------------------------------------------------------------*\
        Class PointParticles Implementation
    \*---------------------------------------------------------------------------*/
    template <typename TypeSystem>
    template <CONCEPT_PARTICLE_CONTAINER STL>
    PointParticles<TypeSystem>::PointParticles(Scalar t, const STL &particle_set) {
        size_t input_num = particle_set.size();
        this->reserve(input_num);
        size_t id = 0;
        for (auto &p : particle_set) {
            pos_.emplace_back(p.pos);
            vel_.emplace_back(p.vel);
            mass_.emplace_back(p.mass);
            idn_.emplace_back(id++);
        }
        time_ = t;
        active_num_ = input_num;
    }

    template <typename TypeSystem>
    void PointParticles<TypeSystem>::resize(size_t new_sz) {
        space::resize_all(new_sz, pos_, vel_, mass_, idn_);
        active_num_ = new_sz;
    }

    template <typename TypeSystem>
    void PointParticles<TypeSystem>::reserve(size_t new_cap) {
        space::reserve_all(new_cap, pos_, vel_, mass_, idn_);
    }

    template <typename TypeSystem>
    void PointParticles<TypeSystem>::clear() {
        space::clear_all(pos_, vel_, mass_, idn_);
        active_num_ = 0;
    }

    template <typename TypeSystem>
    void PointParticles<TypeSystem>::emplace_back(typename PointParticles<TypeSystem>::Particle const &new_particle) {
        pos_.emplace_back(new_particle.pos);
        vel_.emplace_back(new_particle.vel);
        mass_.emplace_back(new_particle.mass);
        idn_.emplace_back(this->number());
        active_num_++;
    }

    template <typename TypeSystem>
    size_t PointParticles<TypeSystem>::number() const {
        return active_num_;
    }

    template <typename TypeSystem>
    size_t PointParticles<TypeSystem>::capacity() const {
        return idn_.capacity();
    }

    template <typename TypeSystem>
    std::string PointParticles<TypeSystem>::column_names() const {
        return "time,id,mass,px,py,pz,vx,vy,vz";
    }

    template <typename TypeSystem>
    std::ostream &operator<<(std::ostream &os, PointParticles<TypeSystem> const &ps) {
        size_t num = ps.number();
        for (size_t i = 0; i < num; ++i) {
            space::print(os, ps.time(), ',', ps.idn(i), ',', ps.mass(i), ',', ps.pos(i), ',', ps.vel(i), '\n');
        }
        return os;
    }
}  // namespace space::particle_set