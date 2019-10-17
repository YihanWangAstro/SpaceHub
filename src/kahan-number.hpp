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
#ifndef SPACEHUB_KAHAN_NUMBER_HPP
#define SPACEHUB_KAHAN_NUMBER_HPP

#include "math.hpp"

namespace space {
/** Kahan number
 *
 *  A way to reduce the round off error when adding a small number to a big one.
 *  See details in [Kahan summation](https://en.wikipedia.org/wiki/Kahan_summation_algorithm)
 */
  template<typename T>
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
    Kahan(T r) : real(r), err(0) {};

    /**
     * Copy constructor.
     * @param k
     */
    Kahan(const Kahan &k) : real(k.real), err(k.err) {};

    /**
     * Assignment operator.
     */
    inline Kahan &operator=(const Kahan &hs) {
      real = hs.real, err = 0;
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
    friend inline Kahan operator-(const Kahan &hs) { return Kahan(-hs.real); }

    /**
     * Addition assignment operator.
     */
    friend inline Kahan &operator+=(Kahan &lhs, const Kahan &rhs) {
      T add = rhs.real - lhs.err;
      T sum = lhs.real + add;

      if (math::abs(add) < math::abs(lhs.real))
        lhs.err = (sum - lhs.real) - add;
      else
        lhs.err = (sum - add) - lhs.real;

      lhs.real = sum;
      return lhs;
    }

    /**
     * Subtraction assignment operator.
     */
    friend inline Kahan &operator-=(Kahan &lhs, const Kahan &rhs) {
      T add = -rhs.real - lhs.err;
      T sum = lhs.real + add;

      if (math::abs(add) < math::abs(lhs.real))
        lhs.err = (sum - lhs.real) - add;
      else
        lhs.err = (sum - add) - lhs.real;

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

  using precise_d = Kahan<double>;
  using precise_f = Kahan<float>;
}  // namespace space
#endif /* kahanNumber_h */
