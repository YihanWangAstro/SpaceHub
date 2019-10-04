#ifndef GENVECT_H
#define GENVECT_H

#include <math.h>
#include <iostream>

namespace space {
/** @brief Self 3D vector class */

template <typename T>
struct Vec3 {
 public:
  /* Typedef */
  using value_type = T;
  /* Typedef */

  value_type x{0};
  value_type y{0};
  value_type z{0};

  Vec3() {}

  Vec3(value_type s) : x(s), y(s), z(s) {}

  Vec3(value_type vx, value_type vy, value_type vz) : x(vx), y(vy), z(vz) {}

  Vec3(const Vec3 &v) : x(v.x), y(v.y), z(v.z) {}

  /** @brief Addition by wise */
  inline Vec3 operator+(const Vec3 &v) const { return Vec3(x + v.x, y + v.y, z + v.z); }

  /** @brief Subtraction by wise */
  inline Vec3 operator-(const Vec3 &v) const { return Vec3(x - v.x, y - v.y, z - v.z); }

  /** @brief Product by wise */
  inline Vec3 operator*(const Vec3 &v) const { return Vec3(x * v.x, y * v.y, z * v.z); }

  /** @brief Divition by wise */
  inline Vec3 operator/(const Vec3 &v) const { return Vec3(x / v.x, y / v.y, z / v.z); }

  /** @brief Add scalar by wise */
  inline Vec3 operator+(const value_type c) const { return Vec3(x + c, y + c, z + c); }

  /** @brief Subtract scalar by wise */
  inline Vec3 operator-(const value_type c) const { return Vec3(x - c, y - c, z - c); }

  /** @brief Multiply scalar by wise */
  inline Vec3 operator*(const value_type c) const { return Vec3(x * c, y * c, z * c); }

  /** @brief Divide scalar by wise */
  inline Vec3 operator/(const value_type c) const { return Vec3(x / c, y / c, z / c); }

  /** @brief Opposite vector */
  inline Vec3 operator-() const { return Vec3(-x, -y, -z); }

  /** @brief Absolute value by wise */
  inline Vec3 abs() const { return Vec3(x > 0 ? x : -x, y > 0 ? y : -y, z > 0 ? z : -z); }

  inline const Vec3 &operator+=(const Vec3 &v) {
    x += v.x, y += v.y, z += v.z;
    return *this;
  }

  inline const Vec3 &operator-=(const Vec3 &v) {
    x -= v.x, y -= v.y, z -= v.z;
    return *this;
  }

  inline const Vec3 &operator*=(const Vec3 &v) {
    x *= v.x, y *= v.y, z *= v.z;
    return *this;
  }

  inline const Vec3 &operator/=(const Vec3 &v) {
    x /= v.x, y /= v.y, z /= v.z;
    return *this;
  }

  inline const Vec3 &operator+=(const value_type c) {
    x += c, y += c, z += c;
    return *this;
  }

  inline const Vec3 &operator-=(const value_type c) {
    x -= c, y -= c, z -= c;
    return *this;
  }

  inline const Vec3 &operator*=(const value_type c) {
    x *= c, y *= c, z *= c;
    return *this;
  }

  inline const Vec3 &operator/=(const value_type c) {
    x /= c, y /= c, z /= c;
    return *this;
  }

  inline Vec3 &operator=(const value_type s) {
    x = y = z = s;
    return *this;
  }

  inline Vec3 &operator=(const Vec3 &v) {
    x = v.x, y = v.y, z = v.z;
    return *this;
  }

  /** @brief Calculate the norm */
  inline value_type norm() const { return sqrt(x * x + y * y + z * z); }

  /** @brief Calculate the norm */
  inline value_type norm2() const { return (x * x + y * y + z * z); }

  inline value_type max_component() const {
    value_type max = (x > y ? x : y);
    return max > z ? max : z;
  }

  /** @brief Calculate the inverse of the norm */
  inline value_type re_norm() const { return 1.0 / sqrt(x * x + y * y + z * z); }

  friend Vec3 operator+(const value_type c, const Vec3 &v) { return Vec3(v.x + c, v.y + c, v.z + c); }

  friend Vec3 operator-(const value_type c, const Vec3 &v) { return Vec3(c - v.x, c - v.y, c - v.z); }

  friend Vec3 operator*(const value_type c, const Vec3 &v) { return Vec3(v.x * c, v.y * c, v.z * c); }

  friend Vec3 operator/(const value_type c, const Vec3 &v) { return Vec3(c / v.x, c / v.y, c / v.z); }

  /** @brief Output to ostream */
  friend std::ostream &operator<<(std::ostream &output, const Vec3 &v) {
    output << v.x << ',' << v.y << ',' << v.z;
    return output;
  }

  /** @brief Input from istream */
  friend std::istream &operator>>(std::istream &input, Vec3 &v) {
    input >> v.x >> v.y >> v.z;
    return input;
  }
};

/** @brief Calculate the Euclid distance of two vectors */
template <typename T>
inline T distance(const Vec3<T> &v1, const Vec3<T> &v2) {
  return sqrt((v1.x - v2.x) * (v1.x - v2.x) + (v1.y - v2.y) * (v1.y - v2.y) + (v1.z - v2.z) * (v1.z - v2.z));
}

/** @brief Calculate the inner product of two vectors */
template <typename T>
inline T dot(const Vec3<T> &v1, const Vec3<T> &v2) {
  return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

template <typename T>
inline T norm(const Vec3<T> &v) {
  return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

template <typename T>
inline T re_norm(const Vec3<T> &v) {
  return 1.0 / sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

template <typename T>
inline T norm2(const Vec3<T> &v) {
  return v.x * v.x + v.y * v.y + v.z * v.z;
}

/** @brief Calculate the cross product of two vectors */
template <typename T>
inline Vec3<T> cross(const Vec3<T> &v1, const Vec3<T> &v2) {
  return Vec3<T>(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
}

using vec3ld = Vec3<long double>;
using vec3d = Vec3<double>;
using vec3f = Vec3<float>;
using vec3i = Vec3<int>;
using vec3c = Vec3<char>;
using vec3b = Vec3<bool>;
}  // namespace space

#include "vector3d.h"   //Specilization of Vec3<double> with AVX;
#include "vector3pd.h"  //Specilization of Vec3<precise_d> with AVX;

#endif