//
// Created by yihan on 4/25/19.
//

#ifndef SPACEHUB_COORDS_HPP
#define SPACEHUB_COORDS_HPP

#include "vector/vector3.h"

namespace space {

/*---------------------------------------------------------------------------*\
    Class Coords Declaration
\*---------------------------------------------------------------------------*/
template <typename T>
struct Coords {
  // type members
  using Scalar = typename T::value_type;

  using Vector = Vec3<Scalar>;

  // constructors
  SPACEHUB_MAKE_CONSTRUCTORS(Coords, default, default, default, default, default);
  
  /**
   *
   * @param count
   */
  explicit Coords(size_t count);
  
  // public methods
  /**
   *
   */
  void clear();

  /**
   *
   * @tparam GenVector
   * @param v
   */
  template <typename GenVector>
  void emplace_back(GenVector const &v);

  /**
   *
   * @param xx
   * @param yy
   * @param zz
   */
  void emplace_back(Scalar xx, Scalar yy, Scalar zz);

  /**
   *
   * @param new_cap
   */
  void reserve(size_t new_cap);

  /**
   *
   * @param new_sz
   */
  void resize(size_t new_sz);

  /**
   *
   */
  void shrink_to_fit();

  /**
   *
   * @return
   */
  size_t size() const;

  // public members
  T x;
  T y;
  T z;
};

/*---------------------------------------------------------------------------*\
    Class Coords Implementation
\*---------------------------------------------------------------------------*/

template <typename T>
Coords<T>::Coords(size_t count) : x(count), y(count), z(count) {}

template <typename T>
size_t Coords<T>::size() const {
  return x.size();
}

template <typename T>
void Coords<T>::reserve(size_t new_cap) {
  reserve_all(new_cap, x, y, z);
}

template <typename T>
void Coords<T>::resize(size_t new_sz) {
  resize_all(new_sz, x, y, z);
}

template <typename T>
template <typename GenVector>
void Coords<T>::emplace_back(GenVector const &v) {
  x.emplace_back(v.x);
  y.emplace_back(v.y);
  z.emplace_back(v.z);
}

template <typename T>
void Coords<T>::emplace_back(Scalar xx, Scalar yy, Scalar zz) {
  x.emplace_back(xx);
  y.emplace_back(yy);
  z.emplace_back(zz);
}

template <typename T>
void Coords<T>::shrink_to_fit() {
  shrink_to_fit_all(x, y, z);
}

template <typename T>
void Coords<T>::clear() {
  clear_all(x, y, z);
}

/*---------------------------------------------------------------------------*\
    Help functions and tools
\*---------------------------------------------------------------------------*/
template <typename STL, typename T>
void add_coords_to(STL &stl, Coords<T> &coords) {
  stl.reserve(coords.size() * 3 + stl.size());
  std::copy(coords.x.begin(), coords.x.end(), std::back_insert_iterator(stl));
  std::copy(coords.y.begin(), coords.y.end(), std::back_insert_iterator(stl));
  std::copy(coords.z.begin(), coords.z.end(), std::back_insert_iterator(stl));
}

template <typename STLIterator, typename T>
void load_to_coords(STLIterator iter_start, STLIterator iter_end, Coords<T> &coords) {
  size_t len = (iter_end - iter_start)/3;
  coords.resize(len);
  auto iter = iter_start;
  for (auto &xx : coords.x) {
    xx = *iter++;
  }
  for (auto &yy : coords.y) {
    yy = *iter++;
  }
  for (auto &zz : coords.z) {
    zz = *iter++;
  }
}

template <typename T>
auto distance(Coords<T> const &c, size_t i, size_t j) {
  auto dx = c.x[i] - c.x[j];
  auto dy = c.y[i] - c.y[j];
  auto dz = c.z[i] - c.z[j];
  return sqrt(dx * dx + dy * dy + dz * dz);
}
}  // namespace space

#endif  // SPACEHUB_COORDS_HPP
