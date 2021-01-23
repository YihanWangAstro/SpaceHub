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
 * @file vector3pd.hpp
 *
 * Header file.
 */
#pragma once

#ifdef __AVX512F__
#pragma message("Using AVX512 on vectorpd")
#include <x86intrin.h>

#include "../kahanNumber.h"

namespace space {
    /** @brief Specilization of vector3d */
    template <>
    struct Vec3<double_p> {
       public:
        /* Typedef */
        using value_type = double_p;
        /* Typedef */

        Vec3<double> __attribute__((aligned(32))) real;
        Vec3<double> __attribute__((aligned(32))) err;

        Vec3() : {};
        Vec3(double vx, double vy, double vz) : mmvalue(_mm512_set_pd(vx, 0, vy, 0, vz, 0, 0, 0.0)){};
        Vec3(__m256d v) : mmvalue(v){};
        Vec3(const Vec3& v) : mmvalue(v.mmvalue){};
        /** @brief Addition by wise */
        inline Vec3 operator+(const Vec3& v) const { return Vec3(_mm512_add_pd(mmvalue, v.mmvalue)); }
        /** @brief Subtraction by wise */
        inline Vec3 operator-(const Vec3& v) const { return Vec3(_mm512_sub_pd(mmvalue, v.mmvalue)); }
        /** @brief Product by wise */
        inline Vec3 operator*(const Vec3& v) const { return Vec3(_mm512_mul_pd(mmvalue, v.mmvalue)); }
        /** @brief Divition by wise */
        inline Vec3 operator/(const Vec3& v) const { return Vec3(_mm512_div_pd(mmvalue, v.mmvalue)); }
        /** @brief Add scalar by wise */
        inline Vec3 operator+(const double_p c) const { return Vec3(_mm512_add_pd(mmvalue, _mm512_set1_pd(c))); }
        /** @brief Subtract scalar by wise */
        inline Vec3 operator-(const double_p c) const { return Vec3(_mm512_sub_pd(mmvalue, _mm512_set1_pd(c))); }
        /** @brief Multiply scalar by wise */
        inline Vec3 operator*(const double_p c) const { return Vec3(_mm512_mul_pd(mmvalue, _mm512_set1_pd(c))); }
        /** @brief Divide scalar by wise */
        inline Vec3 operator/(const double_p c) const { return Vec3(_mm512_div_pd(mmvalue, _mm512_set1_pd(c))); }
        /** @brief Opposite vector */
        inline Vec3 operator-() const { return Vec3(-x, -y, -z); }
        /** @brief Absolute value by wise */
        inline Vec3 abs() const { return Vec3(x > 0 ? x : -x, y > 0 ? y : -y, z > 0 ? z : -z); }

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

        inline const Vec3& operator+=(const Vec3& v) {
            mmvalue = _mm512_add_pd(mmvalue, v.mmvalue);
            return *this;
        }
        inline const Vec3& operator-=(const Vec3& v) {
            mmvalue = _mm512_sub_pd(mmvalue, v.mmvalue);
            return *this;
        }
        inline const Vec3& operator*=(const Vec3& v) {
            mmvalue = _mm512_mul_pd(mmvalue, v.mmvalue);
            return *this;
        }
        inline const Vec3& operator/=(const Vec3& v) {
            mmvalue = _mm512_div_pd(mmvalue, v.mmvalue);
            return *this;
        }
        inline const Vec3& operator+=(const double_p c) {
            mmvalue = _mm512_add_pd(mmvalue, _mm512_set1_pd(c));
            return *this;
        }
        inline const Vec3& operator-=(const double_p c) {
            mmvalue = _mm512_sub_pd(mmvalue, _mm512_set1_pd(c));
            return *this;
        }
        inline const Vec3& operator*=(const double_p c) {
            mmvalue = _mm512_mul_pd(mmvalue, _mm512_set1_pd(c));
            return *this;
        }
        inline const Vec3& operator/=(const double_p c) {
            mmvalue = _mm512_div_pd(mmvalue, _mm512_set1_pd(c));
            return *this;
        }
        inline const Vec3& operator=(const Vec3& v) {
            mmvalue = v.mmvalue;
            return *this;
        }
        /** @brief Calculate the norm */
        inline double_p norm() const {
            Vec3 product = _mm512_mul_pd(mmvalue, mmvalue);
            return sqrt(product.x + product.y + product.z);
        }
        /** @brief Calculate the norm */
        inline double_p norm2() const {
            Vec3 product = _mm512_mul_pd(mmvalue, mmvalue);
            return product.x + product.y + product.z;
        }

        inline double_p max_component() {
            double_p max = (x > y ? x : y);
            return max > z ? max : z;
        }

        /** @brief Calculate the inverse of the norm */
        inline double_p reNorm() const {
            Vec3 product = _mm512_mul_pd(mmvalue, mmvalue);
            return 1.0 / sqrt(product.x + product.y + product.z);
        }
        inline void setZero() { mmvalue = _mm512_setzero_pd(); }
        friend Vec3 operator+(const double_p c, const Vec3& v) {
            return Vec3(_mm512_add_pd(_mm512_set1_pd(c), v.mmvalue));
        }
        friend Vec3 operator-(const double_p c, const Vec3& v) {
            return Vec3(_mm512_sub_pd(_mm512_set1_pd(c), v.mmvalue));
        }
        friend Vec3 operator*(const double_p c, const Vec3& v) {
            return Vec3(_mm512_mul_pd(_mm512_set1_pd(c), v.mmvalue));
        }
        friend Vec3 operator/(const double_p c, const Vec3& v) {
            return Vec3(_mm512_div_pd(_mm512_set1_pd(c), v.mmvalue));
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
    inline double_p distance<double_p>(const Vec3<double_p>& v1, const Vec3<double_p>& v2) {
        __m512d sub = _mm512_sub_pd(v1.mmvalue, v2.mmvalue);
        Vec3<double_p> product = _mm512_mul_pd(sub, sub);
        return sqrt(product.x + product.y + product.z);
    }
    /** @brief Calculate the inner product of two vectors */
    template <>
    inline double_p dot<double_p>(const Vec3<double_p>& v1, const Vec3<double_p>& v2) {
        Vec3<double_p> product = _mm512_mul_pd(v1.mmvalue, v2.mmvalue);
        return product.x + product.y + product.z;
    }
    /** @brief Calculate the cross product of two vectors */
    template <>
    inline Vec3<double_p> cross<double_p>(const Vec3<double_p>& v1, const Vec3<double_p>& v2) {
        return Vec3<double_p>(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
    }

}  // namespace space
#endif