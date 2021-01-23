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
 * @file vector3d.hpp
 *
 * Header file.
 */
#pragma once

#ifdef __AVX__
#pragma message("Using AVX on vector3d")

#include <x86intrin.h>
namespace space {

    /** @brief Specilization of vector3d */
    template <>
    struct alignas(32) Vec3<double> {
       public:
        /* Typedef */
        using value_type = double;
        /* Typedef */

        union alignas(32) {
            struct alignas(32) {
                double x, y, z;
            };
            __m256d mmvalue;
        };
        // __attribute__((aligned(32)));

        Vec3() : mmvalue(_mm256_setzero_pd()){};
        Vec3(double vx, double vy, double vz) : mmvalue(_mm256_set_pd(0.0, vz, vy, vx)){};
        Vec3(double scalar) : mmvalue(_mm256_set1_pd(scalar)){};
        Vec3(__m256d v) : mmvalue(v){};
        Vec3(const Vec3& v) : mmvalue(v.mmvalue){};
        /** @brief Addition by wise */
        inline Vec3 operator+(const Vec3& v) const { return Vec3(_mm256_add_pd(mmvalue, v.mmvalue)); }
        /** @brief Subtraction by wise */
        inline Vec3 operator-(const Vec3& v) const { return Vec3(_mm256_sub_pd(mmvalue, v.mmvalue)); }
        /** @brief Product by wise */
        inline Vec3 operator*(const Vec3& v) const { return Vec3(_mm256_mul_pd(mmvalue, v.mmvalue)); }
        /** @brief Divition by wise */
        inline Vec3 operator/(const Vec3& v) const { return Vec3(_mm256_div_pd(mmvalue, v.mmvalue)); }
        /** @brief Opposite vector */
        inline Vec3 operator-() const {
            // return Vec3(-x, -y, -z);
            return Vec3(-mmvalue);
        }
        /** @brief Absolute value by wise */
        inline Vec3 abs() const {
            // return Vec3(x > 0 ? x : -x, y > 0 ? y : -y, z > 0 ? z : -z);
            return Vec3(_mm256_max_pd(mmvalue, -mmvalue));
        }
        inline const Vec3& operator+=(const Vec3& v) {
            mmvalue = _mm256_add_pd(mmvalue, v.mmvalue);
            return *this;
        }
        inline const Vec3& operator-=(const Vec3& v) {
            mmvalue = _mm256_sub_pd(mmvalue, v.mmvalue);
            return *this;
        }
        inline const Vec3& operator*=(const Vec3& v) {
            mmvalue = _mm256_mul_pd(mmvalue, v.mmvalue);
            return *this;
        }
        inline const Vec3& operator/=(const Vec3& v) {
            mmvalue = _mm256_div_pd(mmvalue, v.mmvalue);
            return *this;
        }
        inline const Vec3& operator=(const Vec3& v) {
            mmvalue = v.mmvalue;
            return *this;
        }

        template <typename U>
        inline Vec3 operator+(const Vec3<U>& v) const {
            return Vec3(x + v.x, y + v.y, z + v.z);
        }

        /** Subtraction by wise */
        template <typename U>
        inline Vec3 operator-(const Vec3<U>& v) const {
            return Vec3(x - v.x, y - v.y, z - v.z);
        }

        /** Product by wise */
        template <typename U>
        inline Vec3 operator*(const Vec3<U>& v) const {
            return Vec3(x * v.x, y * v.y, z * v.z);
        }

        /** Divition by wise */
        template <typename U>
        inline Vec3 operator/(const Vec3<U>& v) const {
            return Vec3(x / v.x, y / v.y, z / v.z);
        }

        /** operator+ for left scalar operation*/
        friend Vec3 operator+(const value_type c, const Vec3& v) {
            return Vec3(_mm256_add_pd(_mm256_set1_pd(c), v.mmvalue));
        }

        /** operator- for left scalar operation*/
        friend Vec3 operator-(const value_type c, const Vec3& v) {
            return Vec3(_mm256_sub_pd(_mm256_set1_pd(c), v.mmvalue));
        }

        /** operator* for left scalar operation*/
        friend Vec3 operator*(const value_type c, const Vec3& v) {
            return Vec3(_mm256_mul_pd(_mm256_set1_pd(c), v.mmvalue));
        }

        /** operator/ for left scalar operation*/
        friend Vec3 operator/(const value_type c, const Vec3& v) {
            return Vec3(_mm256_div_pd(_mm256_set1_pd(c), v.mmvalue));
        }
        /** @brief Output to ostream */
        friend std::ostream& operator<<(std::ostream& output, const Vec3& v) {
            output << v.x << "," << v.y << "," << v.z;
            return output;
        }
        /** @brief Input from istream */
        friend std::istream& operator>>(std::istream& input, Vec3& v) {
            input >> v.x >> v.y >> v.z;
            return input;
        }
    };

    /** @brief Calculate the Euclid distance of two vectors */
    template <>
    inline double distance<double>(const Vec3<double>& v1, const Vec3<double>& v2) {
        __m256d sub = _mm256_sub_pd(v1.mmvalue, v2.mmvalue);
        Vec3<double> product = _mm256_mul_pd(sub, sub);
        return sqrt(product.x + product.y + product.z);
    }
    /** @brief Calculate the inner product of two vectors */
    template <>
    inline double dot<double>(const Vec3<double>& v1, const Vec3<double>& v2) {
        Vec3<double> product = _mm256_mul_pd(v1.mmvalue, v2.mmvalue);
        return product.x + product.y + product.z;
    }
    /** @brief Calculate the cross product of two vectors */
    template <>
    inline Vec3<double> cross<double>(const Vec3<double>& v1, const Vec3<double>& v2) {
        return Vec3<double>(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
    }

}  // namespace space
#endif
