#ifndef TYPECLASS_H
#define TYPECLASS_H

#include "vector/vector3.h"
#include "dev_tools.h"
#include <array>
#include <vector>

namespace SpaceH {

    constexpr size_t DYNAMICAL = 0;

    template<typename T, size_t S>
    struct Array : public std::array<T, S> {};

    template<typename T>
    struct Array<T, DYNAMICAL> : public std::vector<T> {};

    template<typename Real, size_t arraySize>
    struct TypeSystem {
        constexpr static size_t array_size{arraySize};

        using Scalar      = Real;
        using ScalarArray = Array<Scalar, arraySize>;
        using IntArray    = Array<int, arraySize>;
        using IndexArray  = Array<size_t, arraySize>;
        using Vector      = Vec3<Scalar>;
        using VectorArray = Array<Vector, arraySize>;
    };
}//end namespace SpaceH

#endif
