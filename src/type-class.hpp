#ifndef SPACEHUB_TYPE_CLASS_HPP
#define SPACEHUB_TYPE_CLASS_HPP

#include <vector>
#include "coords.hpp"
#include "vector/vector3.h"

namespace space {
/**
  Type system that is used across the Space Hub system. This type class provide all basic type i.e `Scalar`, `Vector`,
  `ScalarArray`, `Coord` and etc,.

  @tparam Real The basic scalar type, i.e `float`, `double`, etc.
  @tparam TContainer The container type, i.e `std::vector`. Use std::vector as default.
 */
  template<typename Real, template<class...> class TContainer = std::vector>
  struct Types {
  public:
    /**
     * @brief Generic container
     *
     * @tparam T variadic template parameters that compatible to various containers.
     */
    template<typename... T>
    using Container = TContainer<T...>;

    /**
     * Floating point like type cross the system
     */
    using Scalar = Real;

    /**
     * 1-d array with value type `Scalar`. Alias of `Container<Scalar>`.
     */
    using ScalarArray = Container<Scalar>;

    /**
     * 1-d array with value type `int`. Alias of `Container<int>`.
     */
    using IntArray = Container<int>;

    /**
     * 1-d array with value type `size_t`. Alias of `Container<size_t>`.
     */
    using IdxArray = Container<size_t>;

    /**
     * 3-d math vector (x, y, z) with `Scalar` type of x, y, z. Use genetic `Vec3` by default, but can be replaced with
     * any other implementation implements interfaces defined in `Vec3`.
     */
    using Vector = Vec3<Scalar>;

    /**
     * 1-d array with value type `Vector`, Alias of `Container<Vector>`.
     */
    using VectorArray = Container<Vector>;

    /**
     * An **Structure of Array** coordinates system (x, y, z) with x, y, z to be a ScalarArray.
     */
    using Coord = Coords<ScalarArray>;
  };
}  // namespace space

#endif
