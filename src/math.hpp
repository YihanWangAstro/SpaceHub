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
#ifndef SPACEHUB_MATH_HPP
#define SPACEHUB_MATH_HPP

#include <cmath>
#include "dev-tools.hpp"

/**
 * @namespace space::math
 * namespace for math
 */
namespace space::math {
/** @brief Self min()*/
template <typename T1, typename T2>
inline T2 min(T1 const &x, T2 const &y) {
  return x > y ? y : x;
}

/** @brief Self max()*/
template <typename T1, typename T2>
inline T2 max(T1 const &x, T2 const &y) {
  return y > x ? y : x;
}

/** @brief Self abs()*/
template <typename T>
inline T abs(const T &x) {
  return x > 0 ? x : -x;
}

template <typename T>
inline T in_range(T low, T x, T high) {
  T tmp = low > x ? low : x;
  return tmp > high ? high : tmp;
}

/**
 *
 * @tparam T
 * @param x
 * @return
 */
template <typename T>
inline constexpr T step(T x) {
  return static_cast<T>(x > 0);
}

/**
 *
 * @tparam T
 * @param x
 * @return
 */
template <typename T>
inline constexpr T sign(T x) {
  return -1 + 2 * static_cast<T>(x > 0);
}

/**
 *
 * @tparam Dtype
 */
template <typename Dtype>
struct epsilon {
  using value_type = typename space::get_value_type<Dtype>::type;
  constexpr static value_type value = std::numeric_limits<value_type>::epsilon();
};

template <typename T>
constexpr T epsilon_v = epsilon<T>::value;

template <typename T>
inline bool iseq(T x, T y) {
  return fabs(x - y) < epsilon_v<T>;
}

/**
 *
 * @tparam Dtype
 */
template <typename Dtype>
struct max_value {
  using value_type = typename space::get_value_type<Dtype>::type;
  constexpr static value_type value = std::numeric_limits<value_type>::max();
};

/**
 *
 * @tparam Dtype
 */
template <typename Dtype>
struct big_value {
  using value_type = typename space::get_value_type<Dtype>::type;
  constexpr static value_type value = 0.1 * std::numeric_limits<value_type>::max();
};

template <typename T>
inline T karmack_fast_inverse_square_root(T x) {
  return 1 / sqrt(x);
}

template <>
inline float karmack_fast_inverse_square_root<float>(float x) {
  float xhalf = 0.5f * x;
  int i = *(int *)&x;
  // i = 0x5f3759df - (i >> 1);
  i = 0x5f375a86 - (i >> 1);
  x = *(float *)&i;
  x = x * (1.5f - xhalf * x * x);
  // x = x*(1.5f - xhalf*x*x);
  return x;
}

template <>
inline double karmack_fast_inverse_square_root<double>(double x) {
  double xhalf = 0.5f * x;
  long long i = *(long long *)&x;
  i = 0x5fe6eb50c7aa19f9 - (i >> 1);
  x = *(double *)&i;
  x = x * (1.5f - xhalf * x * x);
  // x = x*(1.5f - xhalf*x*x);
  return x;
}

/**
 * @brief find root use bisection method
 *
 * @tparam Fun Type Callable ojbect.
 * @param f Callable object/
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
#endif
