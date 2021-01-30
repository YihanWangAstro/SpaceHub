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
    struct Vec3<double_k> {
    private:
        using T = Vec3<double>;
    public:
        using value_type = double_k;
        union{
            struct  {
                double x, y, z;
            };
            T real;
        };
        T err;

        Vec3() : real{}, err{} {};


        Vec3(T r) : real(r), err(0){};


        Vec3(Vec3 const  &k) : real(k.real), err(k.err){};

        /**
         * Assignment operator.
         */
        inline Vec3 &operator=(Vec3 const &hs) {
            real = hs.real, err = hs.err;
            return *this;
        }


        inline operator T() { return real; }


        inline operator T() const { return real; }


        inline void zero_err() { err = 0; }


        friend inline Vec3 operator-(const Vec3 &hs) { return Vec3(-hs.real); }


        friend inline Vec3 &operator+=(Vec3 &lhs, const Vec3 &rhs) {
            T add = rhs.real - lhs.err;
            T sum = lhs.real + add;

            lhs.err = (sum - lhs.real) - add;

            lhs.real = sum;
            return lhs;
        }

        /**
         * Subtraction assignment operator.
         */
        friend inline Vec3 &operator-=(Vec3 &lhs, const Vec3 &rhs) {
            T add = -rhs.real - lhs.err;
            T sum = lhs.real + add;

            lhs.err = (sum - lhs.real) - add;

            lhs.real = sum;
            return lhs;
        }

        /**
         * Division assignment operator.
         */
        friend inline Vec3 &operator/=(Vec3 &lhs, const Vec3 &rhs) {
            lhs.real /= rhs.real;
            return lhs;
        }

        /**
         * Multiple assignment operator.
         */
        friend inline Vec3 &operator*=(Vec3 &lhs, const Vec3 &rhs) {
            lhs.real *= rhs.real;
            return lhs;
        }

        /**
         * Output stream
         */
        friend std::ostream &operator<<(std::ostream &output, const Vec3 &v) {
            output << v.real;
            return output;
        }

        /**
         * Input stream
         */
        friend std::istream &operator>>(std::istream &input, Vec3 &v) {
            input >> v.real;
            v.err = 0;
            return input;
        }

        Vec3 operator+(double c) const {
            return Vec3(real + c);
        }


        Vec3 operator-(double c) const {
            return Vec3(real - c);
        }


        Vec3 operator*(double c) const {
            return Vec3(real * c);
        }


        Vec3 operator/(double c){
            return Vec3(real / c);
        }


        template<typename U, typename V>
        friend Vec3 operator+(Vec3<U> const& v1, Vec3<V> const& v2){
            return Vec3(static_cast<T>(v1)+static_cast<T>(v2));
        }

        template<typename U, typename V>
        friend Vec3 operator-(Vec3<U> const& v1, Vec3<V> const& v2){
            return Vec3(static_cast<T>(v1)-static_cast<T>(v2));
        }

        template<typename U, typename V>
        friend Vec3 operator*(Vec3<U> const& v1, Vec3<V> const& v2){
            return Vec3(static_cast<T>(v1)*static_cast<T>(v2));
        }

        template<typename U, typename V>
        friend Vec3 operator/(Vec3<U> const& v1, Vec3<V> const& v2){
            return Vec3(static_cast<T>(v1)/static_cast<T>(v2));
        }


        friend Vec3 operator+(double c, Vec3 const& v2){
            return Vec3(c+v2.real);
        }


        friend Vec3 operator-(double c, Vec3 const& v2){
            return Vec3(c-v2.real);
        }


        friend Vec3 operator*(double c, Vec3 const& v2){
            return Vec3(c*v2.real);
        }


        friend Vec3 operator/(double c, Vec3 const& v2){
            return Vec3(c/v2.real);
        }

    };

}  // namespace space
