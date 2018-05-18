////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Filename:vector.h                                                                                                   //
//Author:Yihan Wang                                                                                                   //
//                                                                                                                    //
//                                                                                                                    //
//Description:                                                                                                        //
//  Basic 3-d vector class for all vector operations.                                                                 //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef GAR_H
#define GAR_H
#include "chain.h"
#include "../libs.h"
template <typename DataType, size_t N>
class GAR
{
public:
    /////////////////////////////////Type Define////////////////////////////////////
    typedef DataType                   Scalar;
    typedef vec3<Scalar>               Vector;
    typedef std::array<vec3<Scalar>,N> VectorArray;
    typedef std::array<Scalar,N>       ScalarArray;
    typedef std::array<size_t,N>       IndexArray;
    ///////////////////////////////////Interface////////////////////////////////////
    constexpr static size_t size(){return N;}
    constexpr static size_t volume(){return 9*N + 3;}
    
    VectorArray pos;
    VectorArray vel;
    VectorArray auxiVel;
    Scalar      time{0.0};
    Scalar      bindE{0.0};
    Scalar      omega{0.0};
    
    std::array<Scalar, volume()>& array(){ return reinterpret_cast<std::array<Scalar, volume()>&>(pos);}
    void setZero()
    {
        memset(&pos[0], 0, sizeof(Scalar)*volume());
    }
    Scalar getOmega(const ScalarArray& mass)
    {
        return -getPotentialEnergy(mass, pos);
    }
    void initAddiVariable(ScalarArray& mass)
    {
        bindE   = -getTotalEnergy(mass,pos,vel);
        omega   = getOmega(mass);
        auxiVel = vel;
    }
    void toChain(GAR& chainData, IndexArray& index)
    {
        chain::synChain(pos,     chainData.pos,     index);
        chain::synChain(vel,     chainData.vel,     index);
        chain::synChain(auxiVel, chainData.auxiVel, index);
        chainData.time  = time;
        chainData.omega = omega;
        chainData.bindE = bindE;
    }
    void toCartesian(GAR& cartesian, IndexArray& index)
    {
        chain::synCartesian(pos,     cartesian.pos,     index);
        chain::synCartesian(vel,     cartesian.vel,     index);
        chain::synCartesian(auxiVel, cartesian.auxiVel, index);
        cartesian.time  = time;
        cartesian.omega = omega;
        cartesian.bindE = bindE;
    }
    void moveToCentralMassCoords(ScalarArray& mass)
    {
        MoveToCentralMassCoordinate(mass, pos);
        MoveToCentralMassCoordinate(mass, vel);
        MoveToCentralMassCoordinate(mass, auxiVel);
    }
};

#endif

