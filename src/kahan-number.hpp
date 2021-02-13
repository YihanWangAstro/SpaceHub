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
 * @file kahan-number.hpp
 *
 * Header file.
 */
#pragma once

#include "math.hpp"
#include "mpfr.hpp"
namespace space {
    /** Kahan number
     *
     *  A way to reduce the round off error when adding a small number to a big one.
     *  See details in [Kahan summation](https://en.wikipedia.org/wiki/Kahan_summation_algorithm)
     */
    template <typename T>
    struct Kahan {
       public:
        using value_type = T;

        T real, err;

        /**
         * Default constructor
         */
        Kahan() = default;

        /**
         * Single parameter constructor.
         * @param r Scalar
         */
        constexpr Kahan(T r) : real(r), err(0){};

        /**
         * Copy constructor.
         * @param k
         */
        constexpr Kahan(const Kahan &k) : real(k.real), err(k.err){};

        /**
         * Assignment operator.
         */
        // inline Kahan &operator=(const Kahan &hs) {
        //   real = hs.real, err = hs.err;
        //   return *this;
        //}

        /**
         * Conversion operator. Convert Khan number to Scalar i.e `double`, `float`,...
         */
        inline constexpr operator T() { return real; }

        /**
         * Conversion operator. Convert Khan number to Scalar i.e `double`, `float`,...
         */
        inline constexpr operator T() const { return real; }

        /**
         * Set error to 0.
         */
        inline void zero_err() { err = 0; }

        /**
         * Opposite operator.
         */
        friend inline Kahan operator-(const Kahan &hs) { return Kahan(-hs.real); }

        /**
         * Addition assignment operator.
         */
        friend inline Kahan &operator+=(Kahan &lhs, const Kahan &rhs) {
            T add = rhs.real - lhs.err;
            T sum = lhs.real + add;

            lhs.err = (sum - lhs.real) - add;

            lhs.real = sum;
            return lhs;
        }

        /**
         * Subtraction assignment operator.
         */
        friend inline Kahan &operator-=(Kahan &lhs, const Kahan &rhs) {
            T add = -rhs.real - lhs.err;
            T sum = lhs.real + add;

            lhs.err = (sum - lhs.real) - add;

            lhs.real = sum;
            return lhs;
        }

        /**
         * Division assignment operator.
         */
        friend inline Kahan &operator/=(Kahan &lhs, const Kahan &rhs) {
            lhs.real /= rhs.real;
            return lhs;
        }

        /**
         * Multiple assignment operator.
         */
        friend inline Kahan &operator*=(Kahan &lhs, const Kahan &rhs) {
            lhs.real *= rhs.real;
            return lhs;
        }

        /**
         * Output stream
         */
        friend std::ostream &operator<<(std::ostream &output, const Kahan &v) {
            output << v.real;
            return output;
        }

        /**
         * Input stream
         */
        friend std::istream &operator>>(std::istream &input, Kahan &v) {
            input >> v.real;
            v.err = 0;
            return input;
        }
    };

    template <typename T>
    struct Neumaier {
       public:
        using value_type = T;

        T real, err;

        /**
         * Default constructor
         */
        Neumaier() = default;

        /**
         * Single parameter constructor.
         * @param r Scalar
         */
        Neumaier(T r) : real(r), err(0){};

        /**
         * Copy constructor.
         * @param k
         */
        Neumaier(const Neumaier &k) : real(k.real), err(k.err){};

        /**
         * Assignment operator.
         */
        inline Neumaier &operator=(const Neumaier &hs) {
            real = hs.real, err = hs.err;
            return *this;
        }

        /**
         * Conversion operator. Convert Khan number to Scalar i.e `double`, `float`,...
         */
        inline operator T() { return real; }

        /**
         * Conversion operator. Convert Khan number to Scalar i.e `double`, `float`,...
         */
        inline operator T() const { return real; }

        /**
         * Set error to 0.
         */
        inline void zero_err() { err = 0; }

        /**
         * Opposite operator.
         */
        friend inline Neumaier operator-(const Neumaier &hs) { return Neumaier(-hs.real); }

        /**
         * Addition assignment operator.
         */
        friend inline Neumaier &operator+=(Neumaier &lhs, const Neumaier &rhs) {
            T add = rhs.real + lhs.err;
            T sum = lhs.real + add;

            if (math::abs(add) < math::abs(lhs.real))
                lhs.err = (lhs.real - sum) + add;
            else
                lhs.err = (add - sum) + lhs.real;

            lhs.real = sum;
            return lhs;
        }

        /**
         * Subtraction assignment operator.
         */
        friend inline Neumaier &operator-=(Neumaier &lhs, const Neumaier &rhs) {
            T add = -rhs.real + lhs.err;
            T sum = lhs.real + add;

            if (math::abs(add) < math::abs(lhs.real))
                lhs.err = (lhs.real - sum) + add;
            else
                lhs.err = (add - sum) + lhs.real;

            lhs.real = sum;
            return lhs;
        }

        /**
         * Division assignment operator.
         */
        friend inline Neumaier &operator/=(Neumaier &lhs, const Neumaier &rhs) {
            lhs.real /= rhs.real;
            return lhs;
        }

        /**
         * Multiple assignment operator.
         */
        friend inline Neumaier &operator*=(Neumaier &lhs, const Neumaier &rhs) {
            lhs.real *= rhs.real;
            return lhs;
        }

        /**
         * Output stream
         */
        friend std::ostream &operator<<(std::ostream &output, const Neumaier &v) {
            output << v.real;
            return output;
        }

        /**
         * Input stream
         */
        friend std::istream &operator>>(std::istream &input, Neumaier &v) {
            input >> v.real;
            v.err = 0;
            return input;
        }
    };

    template <typename T>
    struct Klein {
       public:
        using value_type = T;

        T real, err, err_of_err;

        /**
         * Default constructor
         */
        Klein() = default;

        /**
         * Single parameter constructor.
         * @param r Scalar
         */
        Klein(T r) : real{r}, err{0}, err_of_err{0} {};

        /**
         * Copy constructor.
         * @param k
         */
        Klein(const Klein &k) : real(k.real), err(k.err), err_of_err{k.err_of_err} {};

        /**
         * Assignment operator.
         */
        inline Klein &operator=(const Klein &hs) {
            real = hs.real, err = hs.err, err_of_err = hs.err_of_err;
            return *this;
        }

        /**
         * Conversion operator. Convert Khan number to Scalar i.e `double`, `float`,...
         */
        inline operator T() { return real; }

        /**
         * Conversion operator. Convert Khan number to Scalar i.e `double`, `float`,...
         */
        inline operator T() const { return real; }

        /**
         * Set error to 0.
         */
        inline void zero_err() {
            err = 0;
            err_of_err = 0;
        }

        /**
         * Opposite operator.
         */
        friend inline Klein operator-(const Klein &hs) { return Klein(-hs.real); }

        /**
         * Addition assignment operator.
         */
        friend inline Klein &operator+=(Klein &lhs, const Klein &rhs) {
            T add = lhs.err + lhs.err_of_err;
            T sum = rhs.real + add;

            if (math::abs(add) < math::abs(rhs.real))
                lhs.err_of_err = (rhs.real - sum) + add;
            else
                lhs.err_of_err = (add - sum) + rhs.real;

            add = sum;
            sum = lhs.real + add;

            if (math::abs(add) < math::abs(lhs.real))
                lhs.err = (lhs.real - sum) + add;
            else
                lhs.err = (add - sum) + lhs.real;

            lhs.real = sum;
            return lhs;
        }

        /**
         * Subtraction assignment operator.
         */
        friend inline Klein &operator-=(Klein &lhs, const Klein &rhs) {
            T add = lhs.err + lhs.err_of_err;
            T sum = -rhs.real + add;

            if (math::abs(add) < math::abs(rhs.real))
                lhs.err_of_err = (-rhs.real - sum) + add;
            else
                lhs.err_of_err = (add - sum) - rhs.real;

            add = sum;
            sum = lhs.real + add;

            if (math::abs(add) < math::abs(lhs.real))
                lhs.err = (lhs.real - sum) + add;
            else
                lhs.err = (add - sum) + lhs.real;

            lhs.real = sum;
            return lhs;
        }

        /**
         * Division assignment operator.
         */
        friend inline Klein &operator/=(Klein &lhs, const Klein &rhs) {
            lhs.real /= rhs.real;
            return lhs;
        }

        /**
         * Multiple assignment operator.
         */
        friend inline Klein &operator*=(Klein &lhs, const Klein &rhs) {
            lhs.real *= rhs.real;
            return lhs;
        }

        /**
         * Output stream
         */
        friend std::ostream &operator<<(std::ostream &output, const Klein &v) {
            output << v.real;
            return output;
        }

        /**
         * Input stream
         */
        friend std::istream &operator>>(std::istream &input, Klein &v) {
            input >> v.real;
            v.err = 0;
            return input;
        }
    };

    // Neumaier<double>;//Kahan<double>;
    using double_k = Kahan<double>;
    using float_k = Kahan<float>;
    using mpreal_k = Kahan<mpfr::mpreal>;

    using double_p = Neumaier<double>;
    using float_p = Neumaier<float>;

    using double_e = Klein<double>;
    using float_e = Klein<float>;
}  // namespace space
