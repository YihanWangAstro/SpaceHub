#ifndef TYPECLASS_H
#define TYPECLASS_H

#include "vector/vector3.h"
#include "coords.h"
#include <vector>

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
        template<typename ...T>
        using Container = TContainer<T...>;

        using Scalar      = Real;
        using ScalarArray = Container<Scalar>;
        using IntArray    = Container<int>;
        using IdxArray    = Container<size_t>;
        using Vector      = Vec3<Scalar>;
        using VectorArray = Container<Vector>;
        using Coord       = Coords<ScalarArray>;
    };
}//end namespace SpaceH

#endif
