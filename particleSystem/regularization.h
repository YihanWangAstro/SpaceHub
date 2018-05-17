////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Filename:logHSystem.h                                                                                               //
//Author:Yihan Wang                                                                                                   //
//                                                                                                                    //
//                                                                                                                    //
//Description:                                                                                                        //
//                                                                                                                    //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef REGULARIZATION_H
#define REGULARIZATION_H
#include "../libs.h"
template<typename DynamicState>
class logH
{
public:
    typedef typename DynamicState::Scalar Scalar;
    constexpr static size_t size(){return DynamicState::size();}
    
    Scalar getPhysicalPosTime(std::array<Scalar, size()>& mass, DynamicState& dyn, Scalar stepSize)
    {
        return stepSize/(dyn.bindE + getKineticEnergy(mass, dyn.vel));
    }
    Scalar getPhysicalVelTime(std::array<Scalar, size()>& mass, DynamicState& dyn, Scalar stepSize)
    {
        return stepSize/(-getPotentialEnergy(mass, dyn.pos));
    }
};

template<typename DynamicState>
class TTL
{
public:
    typedef typename DynamicState::Scalar Scalar;
    constexpr static size_t size(){return DynamicState::size();}
    
    Scalar getPhysicalPosTime(std::array<Scalar, size()>& mass, DynamicState& dyn, Scalar stepSize)
    {
        return stepSize/dyn.omega;
    }
    Scalar getPhysicalVelTime(std::array<Scalar, size()>& mass, DynamicState& dyn, Scalar stepSize)
    {
        return stepSize/dyn.getOmega(mass);
    }
};

template<typename DynamicState>
class NoRegu
{
public:
    typedef typename DynamicState::Scalar Scalar;
    constexpr static size_t size(){return DynamicState::size();}
    
    Scalar getPhysicalPosTime(std::array<Scalar, size()>& mass, DynamicState& dyn, Scalar stepSize)
    {
        return stepSize;
    }
    Scalar getPhysicalVelTime(std::array<Scalar, size()>& mass, DynamicState& dyn, Scalar stepSize)
    {
        return stepSize;
    }
};
#endif

