#ifndef SPACEHUB_TYPE_CLASS_HPP
#define SPACEHUB_TYPE_CLASS_HPP

#include <vector>
#include "coords.hpp"
#include "vector/vector3.h"

namespace space {
/**
 * @brief Type system that is used across the whole spacehub system. This type class provide all basic type i.e
 * 'Scalar', 'Vector', 'ScalarArray', 'Coord' and etc,.
 *
 * @tparam Real The basic scalar type.
 * @tparam TContainer The container type.
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

    using Scalar = Real;
    using ScalarArray = Container<Scalar>;
    using IntArray = Container<int>;
    using IdxArray = Container<size_t>;
    using Vector = Vec3<Scalar>;
    using VectorArray = Container<Vector>;
    using Coord = Coords<ScalarArray>;
  };
}  // namespace space

#endif
