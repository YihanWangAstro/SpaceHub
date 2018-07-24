

#ifndef PROTOTYPE_h
#define PROTOTYPE_h

#include "vector3.h"
#include <array>
#include <vector>

namespace SpaceH
{
    struct EmptyClass{};
    
    constexpr size_t DYNAMICAL = 0;
    
    template<typename T, size_t S>
    struct ArrayWrapper : public std::array<T,S>
    {
        inline void clear(){};
        inline void reserve(size_t new_cap)
        {
            static_assert(new_cap > S, "Fixed size array cannot be allocated!");
        };
        inline void resize(size_t new_size){};
        //inline void
    };
    
    
    
    template<typename T>
    struct ArrayWrapper<T, DYNAMICAL> : public std::vector<T> {};
    
    template<typename T>
    struct get_value_type
    {
    private:
        /*If U has member::value_type, getValueType<T>(0) will match this function. See details on SFINAE. */
        template<typename U>
        static typename U::value_type check(typename U::value_type);
        
        /*If U doesn't have member::value_type, getValueType<T>(0) will match this function. See details on SFINAE. */
        template<typename U>
        static U check(U);
    public:
        using type = decltype(check<T>(0));
    };
    
    template<typename Dtype, size_t Size>
    struct ProtoType
    {
        constexpr static size_t arraySize{Size};
        
        template<typename T, size_t S>
        using Container      = ArrayWrapper<T, S>;
        using Scalar         = Dtype;
        using Vector         = vec3<Scalar>;
        using VectorArray    = Container<Vector, Size>;
        using ScalarArray    = Container<Scalar, Size>;
        using IntArray       = Container<int, Size>;
        using SizeArray      = Container<size_t, Size>;
        using IndexArray     = SizeArray;
        using ScalarBuffer   = Container<Scalar, DYNAMICAL>;
    };
    
    
    enum       PARTICTYPE     { NEUTRONSTAR, STAR, BLACKHOLE, POINT, NONE = 0 };
    /*enum       EVENTTYPE      { TDE, MERGE, ESCAPE, DISRUPTED, UNEVENTFUL, HVS };
    enum class INTEGRATORTYPE { DKDLEAPFROG, KDKLEAPFROG, SYM4, PEFRL, SYM6, SYM8, SYM10 };
    enum class SYSTEMTYPE     { PLAIN, CHAIN };
    enum class ITERTYPE       { BSITER, SEQITER };*/
    enum class NbodyIO {STD, ACTIVE};
    enum class DATASTRUCT     { PLAIN=0, CHAIN };
}//end namespace SpaceH

#endif /* protoType_h */
