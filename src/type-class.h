#ifndef TYPECLASS_H
#define TYPECLASS_H

#include "vector/vector3.h"
#include "dev-tools.h"
#include <array>
#include <vector>

namespace SpaceH {

    template<typename Real, template<class...> class TContainer = std::vector>
    struct TypeSystem {
        template<typename ...T>
        using Container = TContainer<T...>;

        using Scalar   = Real;
        using SArray   = Container<Scalar>;
        using IArray   = Container<int>;
        using IdxArray = Container<size_t>;
        using Vector   = Vec3<Scalar>;
        using VArray   = Container<Vector>;
    };
}//end namespace SpaceH

#endif
