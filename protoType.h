#ifndef PROTOTYPE_h
#define PROTOTYPE_h

#include "vector/vector3.h"
#include "devTools.h"
#include <array>
#include <vector>

namespace SpaceH {
    struct EmptyClass {
    };

    /**
     *
     */
    constexpr size_t DYNAMICAL = 0;

    /**
     *
     * @tparam T
     * @tparam S
     */
    template<typename T, size_t S>
    struct ArrayWrapper : public std::array<T, S> {
        inline void clear() {};
        inline void reserve(size_t new_cap) {};
        inline void resize(size_t new_size) {};
    };

    /**
     *
     * @tparam T
     */
    template<typename T>
    struct ArrayWrapper<T, DYNAMICAL> : public std::vector<T> {
    };

    /**
     *
     * @tparam Dtype
     * @tparam Size
     */
    template<typename Dtype, size_t Size>
    struct ProtoType {
        constexpr static size_t arraySize{Size};

        template<typename T, size_t S>
        using Container      = ArrayWrapper<T, S>;

        using Scalar         = Dtype;
        using Vector         = vec3<Scalar>;
        using VectorArray    = Container<Vector, Size>;
        using ScalarArray    = Container<Scalar, Size>;
        using ScalarBuffer   = Container<Scalar, DYNAMICAL>;
        using IntArray       = Container<int, Size>;
        using SizeArray      = Container<size_t, Size>;
        using IndexArray     = SizeArray;
    };

    enum PARTICTYPE {
        NEUTRONSTAR, STAR, BLACKHOLE, POINT, NONE = 0
    };
    enum class NbodyIO {
        STD, ACTIVE
    };
    enum class DATASTRUCT {
        PLAIN = 0, CHAIN
    };
}//end namespace SpaceH

#endif /* protoType_h */
