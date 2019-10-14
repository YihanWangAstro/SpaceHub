

#ifndef SPACEHUB_KAHAN_NUMBER_HPP
#define SPACEHUB_KAHAN_NUMBER_HPP

#include "own-math.hpp"

namespace space {
/** @brief Kahan number
 *
 *  A way to reduce the round off error when adding a small number to a big one.
 *  See details in https://en.wikipedia.org/wiki/Kahan_summation_algorithm
 */
  template<typename T>
  struct Kahan {
  public:
    using value_type = T;

    T real, err;

    Kahan() = default;

    Kahan(T r) : real(r), err(0) {};

    Kahan(const Kahan &k) : real(k.real), err(k.err) {};

    inline Kahan &operator=(const Kahan &hs) {
      real = hs.real, err = 0;
      return *this;
    }

    inline operator T() { return real; }

    inline operator T() const { return real; }

    inline void zero_err() { err = 0; }

    friend inline Kahan operator-(const Kahan &hs) { return Kahan(-hs.real); }

    friend inline Kahan &operator+=(Kahan &lhs, const Kahan &rhs) {
      T add = rhs.real - lhs.err;
      T sum = lhs.real + add;

      if (space::abs(add) < space::abs(lhs.real))
        lhs.err = (sum - lhs.real) - add;
      else
        lhs.err = (sum - add) - lhs.real;

      lhs.real = sum;
      return lhs;
    }

    friend inline Kahan &operator-=(Kahan &lhs, const Kahan &rhs) {
      T add = -rhs.real - lhs.err;
      T sum = lhs.real + add;

      if (space::abs(add) < space::abs(lhs.real))
        lhs.err = (sum - lhs.real) - add;
      else
        lhs.err = (sum - add) - lhs.real;

      lhs.real = sum;
      return lhs;
    }

    friend inline Kahan &operator/=(Kahan &lhs, const Kahan &rhs) {
      lhs.real /= rhs.real;
      return lhs;
    }

    friend inline Kahan &operator*=(Kahan &lhs, const Kahan &rhs) {
      lhs.real *= rhs.real;
      return lhs;
    }

    /** @brief Output to ostream */
    friend std::ostream &operator<<(std::ostream &output, const Kahan &v) {
      output << v.real;
      return output;
    }

    /** @brief Input from istream */
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
