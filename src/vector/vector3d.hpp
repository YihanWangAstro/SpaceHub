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
 * @file vector3d.hpp
 *
 * Header file.
 */
#pragma once

#ifdef __AVX__
#pragma message("Using AVX on vector3d")

#include <x86intrin.h>

namespace space {

    /** @brief Specialization of vector3d */
    template <>
    struct alignas(32) Vec3<double> {
       public:
        /* Typedef */
        using value_type = double;
        /* Typedef */

        union alignas(32) {
            struct alignas(32) {
                double x{0}, y{0}, z{0};
            };
            __m256d mm_value;
        };
        // __attribute__((aligned(32)));

        Vec3() : mm_value(_mm256_setzero_pd()){};

        Vec3(double vx, double vy, double vz) : mm_value(_mm256_set_pd(0.0, vz, vy, vx)){};

        Vec3(double scalar) : mm_value(_mm256_set1_pd(scalar)){};

        Vec3(__m256d v) : mm_value(v){};

        Vec3(const Vec3 &v) : mm_value(v.mm_value){};

        /** @brief Addition by wise */
        inline Vec3 operator+(const Vec3 &v) const { return Vec3(_mm256_add_pd(mm_value, v.mm_value)); }

        /** @brief Subtraction by wise */
        inline Vec3 operator-(const Vec3 &v) const { return Vec3(_mm256_sub_pd(mm_value, v.mm_value)); }

        /** @brief Product by wise */
        inline Vec3 operator*(const Vec3 &v) const { return Vec3(_mm256_mul_pd(mm_value, v.mm_value)); }

        /** @brief Division by wise */
        inline Vec3 operator/(const Vec3 &v) const { return Vec3(_mm256_div_pd(mm_value, v.mm_value)); }

        /** @brief Opposite vector */
        inline Vec3 operator-() const {
            // return Vec3(-x, -y, -z);
            return Vec3(-mm_value);
        }

        /** @brief Absolute value by wise */
        inline Vec3 abs() const {
            // return Vec3(x > 0 ? x : -x, y > 0 ? y : -y, z > 0 ? z : -z);
            return Vec3(_mm256_max_pd(mm_value, -mm_value));
        }

        inline const Vec3 &operator+=(const Vec3 &v) {
            mm_value = _mm256_add_pd(mm_value, v.mm_value);
            return *this;
        }

        inline const Vec3 &operator-=(const Vec3 &v) {
            mm_value = _mm256_sub_pd(mm_value, v.mm_value);
            return *this;
        }

        inline const Vec3 &operator*=(const Vec3 &v) {
            mm_value = _mm256_mul_pd(mm_value, v.mm_value);
            return *this;
        }

        inline const Vec3 &operator/=(const Vec3 &v) {
            mm_value = _mm256_div_pd(mm_value, v.mm_value);
            return *this;
        }

        /** Addition assignment for vector*/
        template <typename U>
        inline const Vec3 &operator+=(const Vec3<U> &v) {
            x += v.x, y += v.y, z += v.z;
            return *this;
        }

        /** Subtraction assignment for vector*/
        template <typename U>
        inline const Vec3 &operator-=(const Vec3<U> &v) {
            x -= v.x, y -= v.y, z -= v.z;
            return *this;
        }

        /** Multiple assignment for vector*/
        template <typename U>
        inline const Vec3 &operator*=(const Vec3<U> &v) {
            x *= v.x, y *= v.y, z *= v.z;
            return *this;
        }

        /** Division assignment for vector*/
        template <typename U>
        inline const Vec3 &operator/=(const Vec3<U> &v) {
            x /= v.x, y /= v.y, z /= v.z;
            return *this;
        }

        inline Vec3 &operator=(const Vec3 &v) {
            mm_value = v.mm_value;
            return *this;
        }

        template <typename U>
        inline Vec3 operator+(const Vec3<U> &v) const {
            return Vec3(x + v.x, y + v.y, z + v.z);
        }

        /** Subtraction by wise */
        template <typename U>
        inline Vec3 operator-(const Vec3<U> &v) const {
            return Vec3(x - v.x, y - v.y, z - v.z);
        }

        /** Product by wise */
        template <typename U>
        inline Vec3 operator*(const Vec3<U> &v) const {
            return Vec3(x * v.x, y * v.y, z * v.z);
        }

        /** Division by wise */
        template <typename U>
        inline Vec3 operator/(const Vec3<U> &v) const {
            return Vec3(x / v.x, y / v.y, z / v.z);
        }

        /** operator+ for left scalar operation*/
        friend Vec3 operator+(const value_type c, const Vec3 &v) {
            return Vec3(_mm256_add_pd(_mm256_set1_pd(c), v.mm_value));
        }

        /** operator- for left scalar operation*/
        friend Vec3 operator-(const value_type c, const Vec3 &v) {
            return Vec3(_mm256_sub_pd(_mm256_set1_pd(c), v.mm_value));
        }

        /** operator* for left scalar operation*/
        friend Vec3 operator*(const value_type c, const Vec3 &v) {
            return Vec3(_mm256_mul_pd(_mm256_set1_pd(c), v.mm_value));
        }

        /** operator/ for left scalar operation*/
        friend Vec3 operator/(const value_type c, const Vec3 &v) {
            return Vec3(_mm256_div_pd(_mm256_set1_pd(c), v.mm_value));
        }

        /** @brief Output to ostream */
        friend std::ostream &operator<<(std::ostream &output, const Vec3 &v) {
            print_csv(output, v.x, v.y, v.z);
            return output;
        }

        /** @brief Input from istream */
        friend std::istream &operator>>(std::istream &input, Vec3 &v) {
            input >> v.x >> v.y >> v.z;
            return input;
        }
    };

    /** @brief Calculate the Euclid distance of two vectors */
    template <>
    inline double distance<double>(const Vec3<double> &v1, const Vec3<double> &v2) {
        __m256d sub = _mm256_sub_pd(v1.mm_value, v2.mm_value);
        Vec3<double> product = _mm256_mul_pd(sub, sub);
        return sqrt(product.x + product.y + product.z);
    }

    /** @brief Calculate the inner product of two vectors */
    template <>
    inline double dot<double>(const Vec3<double> &v1, const Vec3<double> &v2) {
        Vec3<double> product = _mm256_mul_pd(v1.mm_value, v2.mm_value);
        return product.x + product.y + product.z;
    }

    /** @brief Calculate the cross product of two vectors */
    template <>
    inline Vec3<double> cross<double>(const Vec3<double> &v1, const Vec3<double> &v2) {
        return Vec3<double>(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
    }

}  // namespace space
#endif
