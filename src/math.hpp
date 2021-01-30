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
 * @file math.hpp
 *
 * Header file.
 */
#pragma once

#include <cmath>

#include "dev-tools.hpp"

/**
 * @namespace space::math
 * namespace for math
 */
namespace space::math {
    /**
     * @brief SpaceHub min
     *
     * @tparam T1 Type of the first argument.
     * @tparam T2 Type of the second argument.
     * @param[in] x The first argument.
     * @param[in] y The second argument.
     * @return T2 The min of x, y;
     */
    template <typename T1, typename T2>
    inline constexpr auto min(T1 const &x, T2 const &y) {
        return x > y ? y : x;
    }

    /**
     * @brief SpaceHub max
     *
     * @tparam T1 Type of the first argument.
     * @tparam T2 Type of the second argument.
     * @param[in] x The first argument.
     * @param[in] y The second argument.
     * @return T2 The max of x, y;
     */
    template <typename T1, typename T2>
    inline constexpr auto max(T1 const &x, T2 const &y) {
        return y > x ? y : static_cast<T2>(x);
    }

    /**
     * @brief SpaceHub abs
     *
     * @tparam T Type of the argument.
     * @param[in] x The input argument.
     * @return T The absolute value.
     */
    template <typename T>
    inline constexpr T abs(const T &x) {
        return x >= 0 ? x : -x;
    }

    /**
     * @brief Truncate a number in a range.
     *
     * @tparam T Type of the truncated variable.
     * @param[in] low The lower limit.
     * @param[in] x The variable need to be truncated.
     * @param[in] high The higher limit.
     * @return T The truncated value.
     */
    template <typename T>
    inline constexpr T in_range(T low, T x, T high) {
        T tmp = low > x ? low : x;
        return tmp > high ? high : tmp;
    }

    /**
     * @brief Step function
     *
     * @image html math/step.png width=400px
     *
     * @tparam T Type of the x, y.
     * @param[in] x
     * @return constexpr T y
     */
    template <typename T>
    inline constexpr T step(T x) {
        return static_cast<T>(x > 0);
    }

    /**
     * @brief Sign function
     *
     * @image html math/sign.png width=400px
     *
     * @tparam T Type of x, y
     * @param[in] x
     * @return constexpr T y
     */
    template <typename T>
    inline constexpr T sign(T x) {
        return -1 + 2 * static_cast<T>(x > 0);
    }

    template <typename Dtype>
    struct epsilon {
        using value_type = typename space::raw_type<Dtype>::type;
        constexpr static value_type value = std::numeric_limits<value_type>::epsilon();
    };

    /**
     * @brief Min numerical limit(epsilon) of a floating point type/
     *
     * @tparam T
     */
    template <typename T>
    constexpr T epsilon_v = epsilon<T>::value;

    /**
     * @brief
     *
     * @tparam T
     * @param x
     * @param y
     * @return true
     * @return false
     */
    template <typename T>
    inline bool iseq(T x, T y) {
        return fabs(x - y) < epsilon_v<T>;
    }

    template <typename Dtype>
    struct max_value {
        using value_type = typename space::raw_type<Dtype>::type;
        constexpr static value_type value = std::numeric_limits<value_type>::max();
    };

    /**
     * @brief Max value of a specific type.
     *
     * @tparam T
     */
    template <typename T>
    constexpr T max_value_v = max_value<T>::value;

    template <typename Dtype>
    struct big_value {
        using value_type = typename space::raw_type<Dtype>::type;
        constexpr static value_type value = 0.01 * std::numeric_limits<value_type>::max();
    };

    /**
     * @brief A big value of a specifi type.
     *
     * 0.1 max_value_v
     *
     * @tparam T
     */
    template <typename T>
    constexpr T big_value_v = big_value<T>::value;

    template <typename T>
    inline T karmack_sqrt_inv(T x) {
        return 1 / sqrt(x);
    }

    template <>
    inline float karmack_sqrt_inv<float>(float x) {
        float xhalf = 0.5f * x;
        int i = *(int *) &x;
        // i = 0x5f3759df - (i >> 1);
        i = 0x5f375a86 - (i >> 1);
        x = *(float *) &i;
        x = x * (1.5f - xhalf * x * x);
        // x = x*(1.5f - xhalf*x*x);
        return x;
    }

    template <>
    inline double karmack_sqrt_inv<double>(double x) {
        double xhalf = 0.5f * x;
        long long i = *(long long *) &x;
        i = 0x5fe6eb50c7aa19f9 - (i >> 1);
        x = *(double *) &i;
        x = x * (1.5f - xhalf * x * x);
        // x = x*(1.5f - xhalf*x*x);
        return x;
    }

    /**
     * @brief find root use bisection method
     *
     * @tparam Fun Type of Callable ojbect.
     * @param f Callable object
     * @param low Lower limit of root range
     * @param high Upper limit of root range
     * @return decltype(f(0)) The root.
     */
    template <typename Fun>
    auto root_bisection(Fun f, decltype(f(0)) low, decltype(f(0)) high) -> decltype(f(0)) {
        using Scalar = decltype(f(0));

        for (; fabs((high - low)) > fabs(high) * math::epsilon<Scalar>::value;) {
            Scalar mid = 0.5 * (high + low);
            if (f(mid) > 0)
                high = mid;
            else
                low = mid;
        }
        return 0.5 * (high + low);
    }

    /**
     * @brief Find root with Newton method.
     *
     * @tparam Fun Type of callable object.
     * @param f Callable object.
     * @return decltype(std::declval<Fun>()(0)) The root.
     */
    template <typename Fun>
    decltype(std::declval<Fun>()(0)) root_newton(Fun f) {
        using Scalar = decltype(f(0));
        constexpr size_t max_iter = 1000;
        Scalar x0 = 0;
        Scalar x = 1;
        for (size_t i = 0; i < max_iter; ++i) {
            x = x0 - f(x0);
            if (fabs(x - x0) <= math::epsilon<Scalar>::value) {
                break;
            } else {
                x0 = x;
            }
        }
        return x;
    }
}  // namespace space::math
