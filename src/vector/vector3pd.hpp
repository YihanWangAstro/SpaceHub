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

#include "../kahan-number.hpp"

namespace space {
    template <>
    struct alignas(32) Vec3<double_k> {
       public:
        /* Typedef */
        using value_type = double_k;
        /* Typedef */
        // static_assert(false, "double_k vector");

        double x{0};

        double y{0};

        double z{0};

        double x_err{0};

        double y_err{0};

        double z_err{0};

        // SPACEHUB_MAKE_CONSTRUCTORS(Vec3, default, default, default, default, default);
        Vec3() {}

        explicit Vec3(double s) : x(s), y(s), z(s), x_err(0), y_err(0), z_err(0) {}

        Vec3(double vx, double vy, double vz) : x(vx), y(vy), z(vz), x_err(0), y_err(0), z_err(0) {}

        Vec3(Vec3 const &v) : x(v.x), y(v.y), z(v.z), x_err(v.x_err), y_err(v.y_err), z_err(v.z_err) {}

        template <typename U>
        Vec3(Vec3<U> const &v) : x(v.x), y(v.y), z(v.z), x_err(0), y_err(0), z_err(0) {}

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
        inline Vec3 operator+(const double c) const { return Vec3(x + c, y + c, z + c); }

        /** Subtract scalar by wise */
        inline Vec3 operator-(const double c) const { return Vec3(x - c, y - c, z - c); }

        /** Multiply scalar by wise */
        inline Vec3 operator*(const double c) const { return Vec3(x * c, y * c, z * c); }

        /** Divide scalar by wise */
        inline Vec3 operator/(const double c) const { return Vec3(x / c, y / c, z / c); }

        /** Opposite vector */
        inline Vec3 operator-() const { return Vec3(-x, -y, -z); }

        /** Absolute value by wise */
        inline Vec3 abs() const { return Vec3(x > 0 ? x : -x, y > 0 ? y : -y, z > 0 ? z : -z); }

        /** Addition assignment for vector*/
        template <typename U>
        inline const Vec3 &operator+=(const Vec3<U> &v) {
            double add_x = v.x - x_err;
            double add_y = v.y - y_err;
            double add_z = v.z - z_err;

            double sum_x = x + add_x;
            double sum_y = y + add_y;
            double sum_z = z + add_z;

            x_err = (sum_x - x) - add_x;
            y_err = (sum_y - y) - add_y;
            z_err = (sum_z - z) - add_z;

            x = sum_x, y = sum_y, z = sum_z;
            return *this;
        }

        /** Subtraction assignment for vector*/
        template <typename U>
        inline const Vec3 &operator-=(const Vec3<U> &v) {
            double add_x = -v.x - x_err;
            double add_y = -v.y - y_err;
            double add_z = -v.z - z_err;

            double sum_x = x + add_x;
            double sum_y = y + add_y;
            double sum_z = z + add_z;

            x_err = (sum_x - x) - add_x;
            y_err = (sum_y - y) - add_y;
            z_err = (sum_z - z) - add_z;

            x = sum_x, y = sum_y, z = sum_z;
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
        inline const Vec3 &operator+=(const double c) {
            double add_x = c - x_err;
            double add_y = c - y_err;
            double add_z = c - z_err;

            double sum_x = x + add_x;
            double sum_y = y + add_y;
            double sum_z = z + add_z;

            x_err = (sum_x - x) - add_x;
            y_err = (sum_y - y) - add_y;
            z_err = (sum_z - z) - add_z;

            x = sum_x, y = sum_y, z = sum_z;
            return *this;
        }

        /** Subtraction assignment for scalar*/
        inline const Vec3 &operator-=(const double c) {
            double add_x = -c - x_err;
            double add_y = -c - y_err;
            double add_z = -c - z_err;

            double sum_x = x + add_x;
            double sum_y = y + add_y;
            double sum_z = z + add_z;

            x_err = (sum_x - x) - add_x;
            y_err = (sum_y - y) - add_y;
            z_err = (sum_z - z) - add_z;

            x = sum_x, y = sum_y, z = sum_z;
            return *this;
        }

        /** Multiple assignment for scalar*/
        inline const Vec3 &operator*=(const double c) {
            x *= c, y *= c, z *= c;
            return *this;
        }

        /** Division assignment for scalar*/
        inline const Vec3 &operator/=(const double c) {
            x /= c, y /= c, z /= c;
            return *this;
        }

        /** Assignment operator for scalar*/
        inline Vec3 &operator=(const double s) {
            x = y = z = s;
            x_err = y_err = z_err = 0;
            return *this;
        }

        /** @deprecated Make it non-member function.*/
        inline double norm() const { return sqrt(x * x + y * y + z * z); }

        /** @deprecated Make it non-member function.*/
        inline double norm2() const { return (x * x + y * y + z * z); }

        /** @deprecated Make it non-member function.*/
        inline double max_component() const {
            double max = (x > y ? x : y);
            return max > z ? max : z;
        }

        /** @deprecated Make it non-member function*/
        inline double re_norm() const { return 1.0 / sqrt(x * x + y * y + z * z); }

        /** operator+ for left scalar operation*/
        friend Vec3 operator+(const double c, const Vec3 &v) { return Vec3(v.x + c, v.y + c, v.z + c); }

        /** operator- for left scalar operation*/
        friend Vec3 operator-(const double c, const Vec3 &v) { return Vec3(c - v.x, c - v.y, c - v.z); }

        /** operator* for left scalar operation*/
        friend Vec3 operator*(const double c, const Vec3 &v) { return Vec3(v.x * c, v.y * c, v.z * c); }

        /** operator/ for left scalar operation*/
        friend Vec3 operator/(const double c, const Vec3 &v) { return Vec3(c / v.x, c / v.y, c / v.z); }

        /** output stream */
        friend std::ostream &operator<<(std::ostream &output, const Vec3 &v) {
            output << v.x << ',' << v.y << ',' << v.z;
            return output;
        }

        /** input stream */
        friend std::istream &operator>>(std::istream &input, Vec3 &v) {
            input >> v.x >> v.y >> v.z;
            v.x_err = 0, v.y_err = 0, v.z_err = 0;
            return input;
        }
    };

}  // namespace space
