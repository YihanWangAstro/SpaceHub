
#ifndef LIBS_H
#define LIBS_H

#include "macros.hpp"
#include "own-math.hpp"

namespace space::calc {
template <typename... Args>
constexpr auto add(Args &&... args) {
  return (... + args);
}

template <typename... Args>
constexpr auto mul(Args &&... args) {
  return (... * args);
}

template <typename... Args>
constexpr auto any(Args... args) {
  return (... || args);
}

template <typename... Args>
constexpr auto all(Args... args) {
  return (... && args);
}

template <typename Array>
void array_set_zero(Array &array) {
  for (auto &a : array) {
    a = 0;
  }
}

template <typename... Args>
void set_arrays_zero(Args &... args) {
  (..., (array_set_zero(args)));
}

template <typename Array, typename... Args>
auto array_dot(Array const &a, Array const &b, Args const &... args) {
  DEBUG_MODE_ASSERT(b.size() == a.size(), "length of the array mismatch!");
  typename Array::value_type product{0};
  size_t size = a.size();
  for (size_t i = 0; i < size; ++i) {
    product += (args[i] * ... * (a[i] * b[i]));
  }
  return product;
}

template <typename Array, typename... Args>
void array_add(Array &dst, Array const &a, Array const &b, Args const &... args) {
  DEBUG_MODE_ASSERT(b.size() == a.size() || dst.size() >= a.size(), "length of the array mismatch!");
  size_t size = dst.size();

  for (size_t i = 0; i < size; i++) {
    dst[i] = a[i] + b[i] + (args[i] + ...);
  }
}

template <typename Array, typename... Args>
void array_mul(Array &dst, Array const &a, Array const &b, Args const &... args) {
  size_t size = dst.size();

  for (size_t i = 0; i < size; i++) {
    dst[i] = a[i] * b[i] * (args[i] * ...);
  }
}

template <typename Array>
auto array_sum(Array const &array) {
  typename Array::value_type total = 0;
  for (auto const &a : array) {
    total += a;
  }
  return total;
}

template <typename Scalar, typename Array>
void array_advance(Array &var, Array const &increase, Scalar step_size) {
  size_t size = var.size();

  for (size_t i = 0; i < size; i++) {
    var[i] += increase[i] * step_size;
  }
}

template <typename Array, typename Coord>
void coord_dot(Array &dst, Coord const &a, Coord const &b) {
  size_t size = dst.size();
  for (size_t i = 0; i < size; ++i) {
    dst[i] = a.x[i] * b.x[i] + a.y[i] * b.y[i] + a.z[i] * b.z[i];
  }
}

template <typename Array, typename Coord>
auto coord_contract_to_scalar(Array &coef, Coord const &a, Coord const &b) {
  size_t size = coef.size();
  typename Coord::Scalar sum{0};

  for (size_t i = 0; i < size; ++i) {
    sum += (a.x[i] * b.x[i] + a.y[i] * b.y[i] + a.z[i] * b.z[i]) * coef[i];
  }
  return sum;
}

/*template<typename Scalar, typename Coord>
auto coord_contract_to_scalar(Scalar coef, Coord const &a, Coord const &b) {
    size_t size = a.size();
    Scalar sum{0};

    for (size_t i = 0; i < size; ++i) {
        sum += (a.x[i]*b.x[i] + a.y[i]*b.y[i] + a.z[i]*b.z[i])*coef;
    }
    return sum;
}*/

template <typename Coord>
auto coord_contract_to_scalar(Coord const &a, Coord const &b) {
  size_t size = a.size();
  typename Coord::Scalar sum{0};

  for (size_t i = 0; i < size; ++i) {
    sum += (a.x[i] * b.x[i] + a.y[i] * b.y[i] + a.z[i] * b.z[i]);
  }
  return sum;
}

template <typename Coord, typename... Args>
inline void coord_add(Coord &dst, Coord const &a, Coord const &b, Args const &... args) {
  array_add(dst.x, a.x, b.x, args.x...);
  array_add(dst.y, a.y, b.y, args.y...);
  array_add(dst.z, a.z, b.z, args.z...);
}

template <typename Scalar, typename Coord>
inline void coord_advance(Coord &var, Coord const &increment, Scalar step_size) {
  array_advance(var.x, increment.x, step_size);
  array_advance(var.y, increment.y, step_size);
  array_advance(var.z, increment.z, step_size);
}

template <typename Array>
inline auto calc_com(Array const &mass, Array const &var) {
  return array_dot(var, mass) / array_sum(mass);
}

template <typename Array>
inline auto calc_com(Array const &mass, Array const &var, typename Array::value_type tot_mass) {
  return array_dot(var, mass) / tot_mass;
}

template <typename Array>
void move_to_com(Array &var, typename Array::value_type const &com_var) {
  for (auto &v : var) v -= com_var;
}

template <typename Array1, typename Array2>
inline void move_to_com(Array1 const &mass, Array2 &var) {
  auto com_var = calc_com(mass, var);
  move_to_com(var, com_var);
}

template <typename Coord, typename ScalarArray>
inline void coord_move_to_com(ScalarArray const &mass, Coord &var) {
  auto tot_mass = array_sum(mass);
  move_to_com(var.x, calc_com(mass, var.x, tot_mass));
  move_to_com(var.y, calc_com(mass, var.y, tot_mass));
  move_to_com(var.z, calc_com(mass, var.z, tot_mass));
}

template <typename Particles>
inline auto calc_kinetic_energy(Particles const &ptc) {
  return 0.5 * coord_contract_to_scalar(ptc.mass(), ptc.vel(), ptc.vel());
}

template <typename Particles>
auto calc_potential_energy(Particles const &ptc) {
  typename Particles::Scalar p_eng{0};
  size_t const size = ptc.number();
  auto const &m = ptc.mass();
  auto const &v = ptc.vel();
  auto const &p = ptc.pos();

  for (size_t i = 0; i < size; ++i)
    for (size_t j = i + 1; j < size; ++j) {
      auto dx = p.x[i] - p.x[j];
      auto dy = p.y[i] - p.y[j];
      auto dz = p.z[i] - p.z[j];
      p_eng -= m[i] * m[j] / sqrt(dx * dx + dy * dy + dz * dz);
    }
  return p_eng;
}

template <typename Particles>
inline auto calc_total_energy(Particles const &ptc) {
  return calc_potential_energy(ptc) + calc_kinetic_energy(ptc);
}

/** @brief Calculate the minimal fall free time of two particles
 *
 *  @param  mass mass array of particle.
 *  @param  pos  position array of particle.
 *  @return The minimal fall free time of the two particles
 */
template <typename ScalarArray, typename Coord>
inline auto calc_fall_free_time(ScalarArray const &mass, Coord const &pos) {
  using Scalar = typename ScalarArray::value_type;
  size_t size = mass.size();
  Scalar min_fall_free = max_value<Scalar>::value;

  for (size_t i = 0; i < size; i++) {
    for (size_t j = i + 1; j < size; j++) {
      Scalar dx = pos.x[i] - pos.x[j];
      Scalar dy = pos.y[i] - pos.y[j];
      Scalar dz = pos.z[i] - pos.z[j];
      Scalar r = sqrt(dx * dx + dy * dy + dz * dz);
      Scalar fall_free = pow(r, 1.5) / sqrt(mass[i] + mass[j]);

      if (fall_free < min_fall_free) min_fall_free = fall_free;
    }
  }
  return min_fall_free * space::consts::pi * 0.5 / sqrt(2 * space::consts::G);
}

CREATE_STATIC_MEMBER_CHECK(regu_type);

template <typename Particles>
auto calc_step_scale(Particles const &ptc) {
  if constexpr (HAS_STATIC_MEMBER(Particles, regu_type)) {
    return calc_fall_free_time(ptc.mass(), ptc.pos()) * ptc.omega();
  } else {
    return calc_fall_free_time(ptc.mass(), ptc.pos());
  }
}
}  // namespace space::calc
#endif
