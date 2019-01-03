#ifndef TYPECLASS_H
#define TYPECLASS_H

#include "vector/vector3.h"
#include "dev_tools.h"
#include <array>
#include <vector>

namespace SpaceH {

    template<typename Dtype, size_t Capacity>
    struct TypeSystem {
        constexpr static size_t capacity(){return Capacity;}

        template<typename T, size_t S>
        using Array = std::array<T, S>;

        template<typename T>
        using DynArray = std::vector<T>;

        using Scalar      = Dtype;
        using ScalarArray = Array<Scalar, Capacity>;
        using IntArray    = Array<int, Capacity>;
        using IndexArray  = Array<size_t, Capacity>;
    };

}//end namespace SpaceH

#endif
