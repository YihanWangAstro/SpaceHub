////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Filename:vector.h                                                                                                   //
//Author:Yihan Wang                                                                                                   //
//                                                                                                                    //
//                                                                                                                    //
//Description:                                                                                                        //
//  Basic 3-d vector class for all vector operations.                                                                 //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef DYNAMICSTATE_H
#define DYNAMICSTATE_H
#include "vector3.h"
template <typename DataType, size_t N>
class dynamics
{
public:
    /////////////////////////////////Type Define////////////////////////////////////
    typedef DataType                   Scalar;
    typedef vec3<Scalar>               Vector;
    typedef std::array<vec3<Scalar>,N> VectorArray;
    typedef std::array<Scalar,N>       ScalarArray;
    ///////////////////////////////////Interface////////////////////////////////////
    constexpr static size_t size(){return N;}
    constexpr static size_t volume(){return 6*N + 1;}
    
    VectorArray pos;
    VectorArray vel;
    Scalar      time{0.0};
    
    std::array<Scalar, volume()>& array(){ return reinterpret_cast<std::array<Scalar, volume()>&>(*this);}
    void initAddiVariable(ScalarArray& mass){};
    void setZero()
    {
        memset(&pos[0], 0, sizeof(Scalar)*volume());
    }
};
#endif

