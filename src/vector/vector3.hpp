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
 * @file vector3.hpp
 *
 * Header file.
 */
#pragma once

#include <cmath>
#include <iostream>

#include "../dev-tools.hpp"

namespace space {
    /**
     * Generic 3-d vector (x,y,z) with Scalar x, y, z.
     * */
    template <typename T>
    struct Vec3 {
       public:
        /* Typedef */
        using value_type = T;
        /* Typedef */

        /**
         * x component.
         */
        value_type x{0};
        /**
         * y component.
         */
        value_type y{0};

        /**
         * z component.
         */
        value_type z{0};

        SPACEHUB_MAKE_CONSTRUCTORS(Vec3, default, default, default, default, default);

        /**
         * Create a vector (s,s,s)
         * @param[in] s single scalar.
         */
        explicit Vec3(value_type s) : x(s), y(s), z(s) {}

        /**
         * Create a vector from 3 scalar.
         *
         * @param[in] vx x component.
         * @param[in] vy y component.
         * @param[in] vz z component.
         */
        Vec3(value_type vx, value_type vy, value_type vz) : x(vx), y(vy), z(vz) {}

        /** Addition by wise */
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

        /** Divition by wise */
        template <typename U>
        inline Vec3 operator/(const Vec3<U> &v) const {
            return Vec3(x / v.x, y / v.y, z / v.z);
        }

        /** Add scalar by wise */
        inline Vec3 operator+(const value_type c) const { return Vec3(x + c, y + c, z + c); }

        /** Subtract scalar by wise */
        inline Vec3 operator-(const value_type c) const { return Vec3(x - c, y - c, z - c); }

        /** Multiply scalar by wise */
        inline Vec3 operator*(const value_type c) const { return Vec3(x * c, y * c, z * c); }

        /** Divide scalar by wise */
        inline Vec3 operator/(const value_type c) const { return Vec3(x / c, y / c, z / c); }

        /** Opposite vector */
        inline Vec3 operator-() const { return Vec3(-x, -y, -z); }

        /** Absolute value by wise */
        inline Vec3 abs() const { return Vec3(x > 0 ? x : -x, y > 0 ? y : -y, z > 0 ? z : -z); }

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

        /** Addition assignment for scalar*/
        inline const Vec3 &operator+=(const value_type c) {
            x += c, y += c, z += c;
            return *this;
        }

        /** Subtraction assignment for scalar*/
        inline const Vec3 &operator-=(const value_type c) {
            x -= c, y -= c, z -= c;
            return *this;
        }

        /** Multiple assignment for scalar*/
        inline const Vec3 &operator*=(const value_type c) {
            x *= c, y *= c, z *= c;
            return *this;
        }

        /** Division assignment for scalar*/
        inline const Vec3 &operator/=(const value_type c) {
            x /= c, y /= c, z /= c;
            return *this;
        }

        /** Assignment operator for scalar*/
        inline Vec3 &operator=(const value_type s) {
            x = y = z = s;
            return *this;
        }

        /** @deprecated Make it non-member function.*/
        inline value_type norm() const { return sqrt(x * x + y * y + z * z); }

        /** @deprecated Make it non-member function.*/
        inline value_type norm2() const { return (x * x + y * y + z * z); }

        /** @deprecated Make it non-member function.*/
        inline value_type max_component() const {
            value_type max = (x > y ? x : y);
            return max > z ? max : z;
        }

        /** @deprecated Make it non-member function*/
        inline value_type re_norm() const { return 1.0 / sqrt(x * x + y * y + z * z); }

        /** operator+ for left scalar operation*/
        friend Vec3 operator+(const value_type c, const Vec3 &v) { return Vec3(v.x + c, v.y + c, v.z + c); }

        /** operator- for left scalar operation*/
        friend Vec3 operator-(const value_type c, const Vec3 &v) { return Vec3(c - v.x, c - v.y, c - v.z); }

        /** operator* for left scalar operation*/
        friend Vec3 operator*(const value_type c, const Vec3 &v) { return Vec3(v.x * c, v.y * c, v.z * c); }

        /** operator/ for left scalar operation*/
        friend Vec3 operator/(const value_type c, const Vec3 &v) { return Vec3(c / v.x, c / v.y, c / v.z); }

        /** output stream */
        friend std::ostream &operator<<(std::ostream &output, const Vec3 &v) {
            output << v.x << ',' << v.y << ',' << v.z;
            return output;
        }

        /** input stream */
        friend std::istream &operator>>(std::istream &input, Vec3 &v) {
            input >> v.x >> v.y >> v.z;
            return input;
        }
    };

    /** Calculate the Euclid distance of two vectors */
    template <typename T1, typename T2>
    inline T1 distance(const Vec3<T1> &v1, const Vec3<T2> &v2) {
        auto dx = v1.x - v2.x;
        auto dy = v1.y - v2.y;
        auto dz = v1.z - v2.z;
        return sqrt(dx * dx + dy * dy + dz * dz);
    }

    /** Calculate the inner product of two vectors */
    template <typename T1, typename T2>
    inline T1 dot(const Vec3<T1> &v1, const Vec3<T2> &v2) {
        return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
    }

    /** Calculate the length of a vector*/
    template <typename T>
    inline T norm(const Vec3<T> &v) {
        return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    }

    /**  Calculate the inverse lenght of a vector*/
    template <typename T>
    inline T re_norm(const Vec3<T> &v) {
        return 1.0 / sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    }

    /** Calculate the length square of a vector*/
    template <typename T>
    inline T norm2(const Vec3<T> &v) {
        return v.x * v.x + v.y * v.y + v.z * v.z;
    }

    /** Calculate the cross product of two vectors */
    template <typename T1, typename T2>
    inline Vec3<T1> cross(const Vec3<T1> &v1, const Vec3<T2> &v2) {
        return Vec3<T1>(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
    }

    template <typename T>
    inline Vec3<T> vec_abs(Vec3<T> const &v) {
        return Vec3<T>(fabs(v.x), fabs(v.y), fabs(v.z));
    }

    template <typename T>
    inline Vec3<T> vec_max(Vec3<T> const &v1, Vec3<T> const &v2) {
        return Vec3<T>(std::max(v1.x, v2.x), std::max(v1.y, v2.y), std::max(v1.z, v2.z));
    }

    template <typename T>
    inline Vec3<T> vec_min(Vec3<T> const &v1, Vec3<T> const &v2) {
        return Vec3<T>(std::min(v1.x, v2.x), std::min(v1.y, v2.y), std::min(v1.z, v2.z));
    }

    template <typename T>
    inline T max_abs(Vec3<T> const &v) {
        auto max = std::max(fabs(v.x), fabs(v.y));
        return std::max(max, fabs(v.z));
    }

    /**
     * Alias for long double 3-d vector
     */
    using vec3ld = Vec3<long double>;
    /**
     * Alias for double 3-d vector
     */
    using vec3d = Vec3<double>;
    /**
     * Alias for float 3-d vector
     */
    using vec3f = Vec3<float>;
    /**
     * Alias for int 3-d vector
     */
    using vec3i = Vec3<int>;
    /**
     * Alias for char 3-d vector
     */
    using vec3c = Vec3<char>;
    /**
     * Alias for bool 3-d vector
     */
    using vec3b = Vec3<bool>;
}  // namespace space

//#include "vector3d.hpp"  //Specilization of Vec3<double> with AVX;
//#include "vector3pd.hpp"  //Specilization of Vec3<precise_d> with AVX;
