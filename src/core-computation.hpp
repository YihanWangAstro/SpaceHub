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
 * @file core-computation.hpp
 *
 * Header file.
 */
#pragma once
#ifdef _OPENMP
#include <omp.h>
#endif
#include <algorithm>
#include <numeric>
#include <type_traits>

#include "macros.hpp"
#include "math.hpp"
#include "spacehub-concepts.hpp"
#include "vector/vector3.hpp"

#ifdef __AVX__

#include <x86intrin.h>

#endif
/**
 * @namespace space::calc
 * Documentation for space
 */
namespace space::calc {

    /**
     * @brief Calculate the sum of variables
     *
     * @tparam Args The types of variables
     * @param[in] args Variables.
     * @return constexpr auto The sum.
     */
    template <typename... Args>
    constexpr auto add(Args &&...args) {
        return (... + args);
    }

    /**
     * @brief Calculate the product of variables.
     *
     * @tparam Args The types of variables.
     * @param[in] args Variables.
     * @return constexpr auto The product.
     */
    template <typename... Args>
    constexpr auto mul(Args &&...args) {
        return (... * args);
    }

    /**
     * @brief Calculate the *logical or* of variables.
     *
     * @tparam Args The types of variables.
     * @param[in] args Variables.
     * @return constexpr auto The result.
     */
    template <typename... Args>
    constexpr auto any(Args... args) {
        return (... || args);
    }

    /**
     * @brief Calculate the *logical and* of variables.
     *
     * @tparam Args The types of variables.
     * @param[in] args Variables.
     * @return constexpr auto The results.
     */
    template <typename... Args>
    constexpr auto all(Args... args) {
        return (... && args);
    }

    /**
     * @brief Set all elements of a std::ranges(Container) to 0.
     *
     * @tparam Array Type of container.
     * @param[in,out] array Container.
     */
    template <typename Array>
    void array_set_zero(Array &array) {
        std::fill(array.begin(), array.end(), typename Array::value_type{0});
        /* for (auto &a : array) {
             a = 0;
         }*/
    }

    /**
     * @brief Set all elements of std::ranges(Containers) to 0.
     *
     * @tparam Args The types of containers
     * @param[in,out] args Containers.
     */
    template <typename... Args>
    void set_arrays_zero(Args &...args) {
        (..., (array_set_zero(args)));
    }

    /**
     * @brief Calculate the (general)dot product of arrays.
     *
     *  @f$ \Sigma_{i=0}^n a_ib_ic_i...@f$
     *
     * @tparam Array The type of the first two arrays.
     * @tparam Args The types of the arrays.
     * @param[in] a The first array.
     * @param[in] b The second array.
     * @param[in] args The rest arrays.
     * @return auto The dot product.
     *
     * @note One should ensure the input arrays have method `size()` and they have the same length.
     */
    template <typename Array, typename Array2, typename... Args>
    auto array_dot(Array const &a, Array2 const &b, Args const &...args) {
        DEBUG_MODE_ASSERT(b.size() == a.size(), "length of the array mismatch!");
        typename Array::value_type product{0};
        size_t const size = a.size();
        for (size_t i = 0; i < size; ++i) {
            product += (args[i] * ... * (a[i] * b[i]));
        }
        return product;
    }

    /**
     * @brief Element wise addition of arrays.
     *
     * @f$ dst[i] = a[i] + b[i] + c[i] + ...@f$
     *
     * @tparam Array The type of the first array.
     * @tparam Args The type of the arrays.
     * @param[out] dst The destination array to store the result.
     * @param[in] a The first input array.
     * @param[in] args The rest arrays.
     */
    template <typename Array, typename... Args>
    void array_add(Array &dst, Args const &...args) {
        size_t const size = dst.size();

//#pragma omp parallel for
#pragma GCC ivdep
        for (size_t i = 0; i < size; i++) {
            dst[i] = (args[i] + ...);
        }
    }

    CREATE_MEMBER_CHECK(err);
    CREATE_MEMBER_CHECK(x);
    CREATE_MEMBER_CHECK(y);
    CREATE_MEMBER_CHECK(z);

#define IS_VECTOR3(CLASS) HAS_MEMBER(CLASS, x) && HAS_MEMBER(CLASS, y) && HAS_MEMBER(CLASS, z)

    template <typename A1, typename A2>
    void array_load_err(A1 &dst, A2 const &src) {
        if constexpr (HAS_MEMBER(typename A1::value_type, err)) {
            size_t size = dst.size();
            for (size_t i = 0; i < size; ++i) {
                dst[i].err = src[i];
            }
        } else if constexpr (IS_VECTOR3(typename A1::value_type) &&
                             HAS_MEMBER(typename A1::value_type::value_type, err)) {
            size_t size = dst.size();
            for (size_t i = 0; i < size; ++i) {
                dst[i].x.err = src[i].x;
                dst[i].y.err = src[i].y;
                dst[i].z.err = src[i].z;
            }
        }
    }

    template <typename A1, typename A2>
    void array_save_err(A1 &dst, A2 const &src) {
        if constexpr (HAS_MEMBER(typename A2::value_type, err)) {
            size_t size = dst.size();
            for (size_t i = 0; i < size; ++i) {
                dst[i] = src[i].err;
            }
        } else if constexpr (IS_VECTOR3(typename A2::value_type) &&
                             HAS_MEMBER(typename A2::value_type::value_type, err)) {
            size_t size = dst.size();
            for (size_t i = 0; i < size; ++i) {
                dst[i].x = src[i].x.err;
                dst[i].y = src[i].y.err;
                dst[i].z = src[i].z.err;
            }
        }
    }

    template <typename A1, typename A2>
    void array_copy_err(A1 &dst, A2 const &src) {
        if constexpr (HAS_MEMBER(typename A2::value_type, err)) {
            size_t size = dst.size();
            for (size_t i = 0; i < size; ++i) {
                dst[i].err = src[i].err;
            }
        } else if constexpr (IS_VECTOR3(typename A2::value_type) &&
                             HAS_MEMBER(typename A2::value_type::value_type, err)) {
            size_t size = dst.size();
            for (size_t i = 0; i < size; ++i) {
                dst[i].x.err = src[i].x.err;
                dst[i].y.err = src[i].y.err;
                dst[i].z.err = src[i].z.err;
            }
        }
    }

    template <typename A1>
    void array_set_err_zero(A1 &dst) {
        if constexpr (HAS_MEMBER(typename A1::value_type, err)) {
            size_t size = dst.size();
            for (size_t i = 0; i < size; ++i) {
                dst[i].err = 0;
            }
        } else if constexpr (IS_VECTOR3(typename A1::value_type) &&
                             HAS_MEMBER(typename A1::value_type::value_type, err)) {
            size_t size = dst.size();
            for (size_t i = 0; i < size; ++i) {
                dst[i].x.err = 0;
                dst[i].y.err = 0;
                dst[i].z.err = 0;
            }
        }
    }

    /**
     * @brief Element wise scale.
     *
     * @f$ dst[i] = scale * a[i]@f$
     *
     * @tparam Array The type of the array.
     * @tparam Scalar The type of the scale factor.
     * @param[out] dst The destination array to store the result.
     * @param[in] a The input array.
     * @param[in] scale Scale factor.
     */
    template <typename Array1, typename Array2, typename Scalar>
    void array_scale(Array1 &dst, Array2 const &a, Scalar scale) {
        DEBUG_MODE_ASSERT(b.size() == a.size() || dst.size() >= a.size(), "length of the array mismatch!");

        // std::transform(a.begin(), a.end(), dst.begin(), [=](auto x) { return scale * x; });
        size_t const size = dst.size();
//#pragma omp parallel for
#pragma GCC ivdep
        for (size_t i = 0; i < size; i++) {
            dst[i] = a[i] * scale;
        }
    }

    template <typename Array1, typename Array2, typename Array3, typename Scalar>
    void array_scale_add(Array1 &dst, Array2 const &a, Array3 const &b, Scalar c) {
        std::transform(a.begin(), a.end(), b.begin(), dst.begin(), [=](auto x, auto y) { return (x + y) * c; });
    }

    template <typename Array1, typename Array2, typename Array3, typename Scalar>
    void array_scale_add(Array1 &dst, Array2 const &a, Scalar b, Array3 const &c, Scalar d) {
        std::transform(a.begin(), a.end(), c.begin(), dst.begin(), [=](auto x, auto y) { return x * b + y * d; });
    }

    template <typename Array1, typename Array2, typename Scalar>
    void array_div_scale(Array1 &dst, Array2 const &a, Scalar scale) {
        DEBUG_MODE_ASSERT(b.size() == a.size() || dst.size() >= a.size(), "length of the array mismatch!");
        // std::transform(a.begin(), a.end(), dst.begin(), [=](auto x) { return x / scale; });
        size_t const size = dst.size();
//#pragma omp parallel for
#pragma GCC ivdep
        for (size_t i = 0; i < size; i++) {
            dst[i] = a[i] / scale;
        }
    }

    /**
     * @brief Element wise multiplication of arrays.
     *
     * @f$ dst[i] = a[i] * b[i] * c[i] * ...@f$
     *
     * @tparam Array The type of the first array.
     * @tparam Args The type of the arrays.
     * @param[out] dst The destination array to store the result.
     * @param[in] a The first input array.
     * @param[in] b The second input array.
     * @param[in] args The rest arrays.
     */
    template <typename Array, typename... Args>
    void array_mul(Array &dst, Args const &...args) {
        size_t const size = dst.size();
//#pragma omp parallel for
#pragma GCC ivdep
        for (size_t i = 0; i < size; i++) {
            dst[i] = (args[i] * ...);
        }
    }

    template <typename Array1, typename Array2, typename Array3>
    void array_sub(Array1 &dst, Array2 const &a, Array3 const &b) {
        // DEBUG_MODE_ASSERT(b.size() == a.size() || dst.size() >= a.size(), "length of the array mismatch!");
        // std::transform(a.begin(), a.end(), b.begin(), dst.begin(), [](auto x, auto y) { return x - y; });
        size_t const size = dst.size();
//#pragma omp parallel for
#pragma GCC ivdep
        for (size_t i = 0; i < size; i++) {
            dst[i] = a[i] - b[i];
        }
    }

    template <typename Array1, typename Array2, typename Array3, typename Scalar>
    void array_scale_sub(Array1 &dst, Array2 const &a, Array3 const &b, Scalar c) {
        std::transform(a.begin(), a.end(), b.begin(), dst.begin(), [=](auto x, auto y) { return (x - y) * c; });
    }

    template <typename Array1, typename Array2, typename Array3>
    void array_div(Array1 &dst, Array2 const &a, Array3 const &b) {
        // DEBUG_MODE_ASSERT(b.size() == a.size() || dst.size() >= a.size(), "length of the array mismatch!");
        size_t const size = dst.size();
//#pragma omp parallel for
#pragma GCC ivdep
        for (size_t i = 0; i < size; i++) {
            dst[i] = a[i] / b[i];
        }
    }

    template <typename Array>
    auto array_sum(Array const &array) {
        // return std::reduce(array.begin(), array.end());
        typename Array::value_type total = 0;

        for (auto const &a : array) {
            total += a;
        }
        return total;
    }

    template <typename Array1, typename Array2>
    void array_advance(Array1 &var, Array2 const &increment) {
        size_t const size = var.size();
#pragma GCC ivdep
        for (size_t i = 0; i < size; i++) {
            var[i] += increment[i];
        }
    }

    template <typename Scalar, typename Array1, typename Array2>
    void array_advance(Array1 &var, Array2 const &increment, Scalar step_size) {
        size_t const size = var.size();
#pragma GCC ivdep
        for (size_t i = 0; i < size; i++) {
            var[i] += increment[i] * step_size;
        }
    }

    template <typename Scalar, typename Array1, typename Array2, typename Array3>
    void array_advance(Array1 &dst, Array2 const &var, Array3 const &increment, Scalar step_size) {
        size_t const size = var.size();
//#pragma omp parallel for
#pragma GCC ivdep
        for (size_t i = 0; i < size; i++) {
            dst[i] = var[i] + increment[i] * step_size;
        }
    }

    template <typename Array1, typename Array2>
    void array_retreat(Array1 &var, Array2 const &increment) {
        size_t const size = var.size();
#pragma GCC ivdep
        for (size_t i = 0; i < size; i++) {
            var[i] -= increment[i];
        }
    }

    template <typename Scalar, typename Array1, typename Array2>
    void array_retreat(Array1 &var, Array2 const &increment, Scalar step_size) {
        size_t const size = var.size();
#pragma GCC ivdep
        for (size_t i = 0; i < size; i++) {
            var[i] -= increment[i] * step_size;
        }
    }

    template <typename Scalar, typename Array1, typename Array2, typename Array3>
    void array_retreat(Array1 &dst, Array2 const &var, Array3 const &increment, Scalar step_size) {
        size_t const size = var.size();
//#pragma omp parallel for
#pragma GCC ivdep
        for (size_t i = 0; i < size; i++) {
            dst[i] = var[i] - increment[i] * step_size;
        }
    }

    template <typename Array, typename VectorArray1, typename VectorArray2>
    void coord_dot(Array &dst, VectorArray1 const &a, VectorArray2 const &b) {
        size_t const size = dst.size();
//#pragma omp parallel for
#pragma GCC ivdep
        for (size_t i = 0; i < size; ++i) {
            dst[i] = dot(a[i], b[i]);
        }
    }

    template <typename Array, typename VectorArray1, typename VectorArray2>
    auto coord_contract_to_scalar(Array &coef, VectorArray1 const &a, VectorArray2 const &b) {
        size_t const size = coef.size();
        typename Array::value_type sum{0};

        for (size_t i = 0; i < size; ++i) {
            sum += dot(a[i], b[i]) * coef[i];
        }
        return sum;
    }

    /*template<typename Scalar, typename Coord>
    auto coord_contract_to_scalar(Scalar table_coef, Coord const &a, Coord const &b) {
        size_t size = a.size();
        Scalar sum{0};

        for (size_t i = 0; i < size; ++i) {
            sum += (a.x[i]*b.x[i] + a.y[i]*b.y[i] + a.z[i]*b.z[i])*table_coef;
        }
        return sum;
    }*/

    template <typename VectorArray1, typename VectorArray2>
    auto coord_contract_to_scalar(VectorArray1 const &a, VectorArray2 const &b) {
        size_t const size = a.size();
        typename VectorArray1::value_type::value_type sum{0};

        for (size_t i = 0; i < size; ++i) {
            sum += dot(a[i], b[i]);
        }
        return sum;
    }

    template <typename ScalarArray, typename VectorArray>
    inline auto calc_com(ScalarArray const &mass, VectorArray const &var) -> typename VectorArray::value_type {
        // typename VectorArray::value_type com = array_dot(var, mass) / array_sum(mass);//evaluate for lazy vecor
        // otherwise will return an expression return com;
        return array_dot(var, mass) / array_sum(mass);
    }

    template <typename Array1, typename Array2>
    inline void move_to_com(Array1 const &mass, Array2 &var) {
        auto com_var = calc_com(mass, var);
#pragma GCC ivdep
        for (auto &v : var) v -= com_var;
    }

    template <CONCEPT_PARTICLES_DATA Particles>
    inline auto calc_kinetic_energy(Particles const &ptc) -> typename Particles::Scalar {
        return 0.5 * coord_contract_to_scalar(ptc.mass(), ptc.vel(), ptc.vel());
    }

    CREATE_METHOD_CHECK(chain_pos);

    CREATE_METHOD_CHECK(index);

    CREATE_METHOD_CHECK(omega);

    CREATE_METHOD_CHECK(bindE);

    CREATE_METHOD_CHECK(mass);

    CREATE_METHOD_CHECK(pos);

    CREATE_METHOD_CHECK(vel);

    CREATE_STATIC_MEMBER_CHECK(regu_type);

    template <CONCEPT_PARTICLES_DATA Particle>
    auto calc_potential_energy(Particle const &particle1, Particle const &particle2) -> typename Particle::Scalar {
        typename Particle::Scalar potential_eng = -consts::G * particle1.mass * particle2.mass;
        return potential_eng / norm(particle1.pos - particle2.pos);
    }

    template <CONCEPT_PARTICLES_DATA Particles>
    auto calc_potential_energy(Particles const &particles) -> typename Particles::Scalar {
        typename Particles::Scalar potential_eng{0};
        size_t const size = particles.number();
        auto const &m = particles.mass();
        auto const &v = particles.vel();
        auto const &p = particles.pos();

        if constexpr (HAS_METHOD(Particles, chain_pos) && HAS_METHOD(Particles, index)) {
            auto const &ch_p = particles.chain_pos();
            auto const &idx = particles.index();

            for (size_t i = 0; i < size - 1; ++i) {
                potential_eng -= m[idx[i]] * m[idx[i + 1]] / norm(ch_p[i]);
            }

            for (size_t i = 0; i < size - 2; ++i) {
                decltype(ch_p[i]) dr = ch_p[i] + ch_p[i + 1];
                potential_eng -= m[idx[i]] * m[idx[i + 2]] / norm(dr);
            }

            for (size_t i = 0; i < size; ++i) {
                for (size_t j = i + 3; j < size; ++j) {
                    decltype(p[idx[j]]) dr = p[idx[j]] - p[idx[i]];
                    potential_eng -= m[idx[i]] * m[idx[j]] / norm(dr);
                }
            }
        } else {
            for (size_t i = 0; i < size; ++i)
                for (size_t j = i + 1; j < size; ++j) {
                    decltype(p[i]) dr = p[i] - p[j];
                    potential_eng -= m[i] * m[j] / norm(dr);
                }
        }
        return potential_eng * consts::G;
    }

    template <CONCEPT_PARTICLES_DATA Particles>
    inline auto calc_total_energy(Particles const &particles) -> typename Particles::Scalar {
        return calc_potential_energy(particles) + calc_kinetic_energy(particles);
    }

    /** @brief Calculate the minimal fall free time of two particles
     *
     *  @param  mass mass array of particle.
     *  @param  position  position array of particle.
     *  @return The minimal fall free time of the two particles
     */
    template <typename ScalarArray, typename VectorArray>
    inline auto calc_fall_free_time(ScalarArray const &mass, VectorArray const &position) {
        using Scalar = typename ScalarArray::value_type;
        size_t const size = mass.size();
        Scalar min_fall_free = math::max_value<double>::value;  // TODO:

        for (size_t i = 0; i < size; i++) {
            for (size_t j = i + 1; j < size; j++) {
                Scalar r = distance(position[i], position[j]);
                Scalar fall_free = POW(r, 1.5) / sqrt(mass[i] + mass[j]);
                if (fall_free < min_fall_free) min_fall_free = fall_free;
            }
        }
        return min_fall_free * space::consts::pi * 0.5 / sqrt(2 * space::consts::G);
    }

    /**
     *
     * @tparam Particles
     * @param particles
     * @return
     */
    template <CONCEPT_PARTICLES_DATA Particles>
    auto calc_step_scale(Particles const &particles) {
        if constexpr (HAS_STATIC_MEMBER(Particles, regu_type)) {
            return particles.step_scale();
        } else {
            return 1.0;
        }
    }

    template <CONCEPT_PARTICLES_DATA Particles>
    auto calc_energy_error(Particles const &particles, typename Particles::Scalar E0) -> typename Particles::Scalar {
        using Scalar = typename Particles::Scalar;
        Scalar U = -calc_potential_energy(particles);
        Scalar T = calc_kinetic_energy(particles);
        if constexpr (HAS_METHOD(Particles, bindE) && HAS_STATIC_MEMBER(Particles, regu_type)) {
            return LOG(fabs((T + particles.bindE()) / U));
        } else {
            return fabs((T - U - E0) / E0);
        }
    }
}  // namespace space::calc
