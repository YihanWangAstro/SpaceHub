////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Filename:libs.h                                                                                                     //
//Author:Yihan Wang                                                                                                   //
//                                                                                                                    //
//                                                                                                                    //
//Description:                                                                                                        //
//                                                                                                                    //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef LIBS_H
#define LIBS_H
#include "vector3.h"
#include "rdFloat.h"

template<typename T1, typename T2>
inline const T2 min(const T1& x, const T2& y)
{
    return x > y ? y : x;
}

template<typename T1, typename T2>
inline const T2 max(const T1& x, const T2& y)
{
    return y > x ? y : x;
}

template<class T>
inline const T abs(const T& x)
{
    return x > -x ? x : -x;
}

template <class T>
void swap(T& a, T& b)
{
    T tmp(std::move(a));
    a = std::move(b);
    b = std::move(tmp);
}

template<typename Scalar, size_t N>
void KahanAdvance(std::array<vec3<Scalar>, N>& var, const std::array<vec3<Scalar>, N>& increase, std::array<vec3<Scalar>, N>& err, Scalar dt)
{
    vec3<Scalar> add;
    vec3<Scalar> sum;
    for(size_t i = 0 ; i < N; ++i)
    {
        add    = increase[i]*dt - err[i];
        sum    = var[i] + add;
        err[i] = (sum - var[i]) - add;
        var[i] = sum;
    }
}

template<typename Scalar>
void KahanAdvance(Scalar& var, const Scalar increase, Scalar& err)
{
    Scalar add;
    Scalar sum;

    add = increase - err;
    sum = var + add;
    err = (sum - var) - add;
    var = sum;
}

template<typename Scalar, size_t N>
void advanceVariable(std::array<vec3<Scalar>, N>& var, const std::array<vec3<Scalar>, N>& add, Scalar dt)
{
    for(size_t i = 0 ; i < N; ++i)
        var[i] += add[i]*dt;
}

template<typename Scalar, size_t N>
void MoveToCentralMassCoordinate(const std::array<Scalar,N>& mass, std::array<vec3<Scalar>, N>& phyVar)
{
    vec3<Scalar> centralMassVar(0.0, 0.0, 0.0);
    Scalar totalMass = 0;
    
    for(size_t i = 0 ; i < N ; ++i)
    {
        totalMass += mass[i];
        centralMassVar += phyVar[i]*mass[i];
    }
    
    centralMassVar /= totalMass;
   
    for(size_t i = 0 ; i < N ; ++i)
    {
        phyVar[i] -= centralMassVar;
    }
}

template<typename Scalar, size_t N>
double getKineticEnergy(const std::array<Scalar,N>& mass, const std::array<vec3<Scalar>, N>& vel)
{
    Scalar kineticEnergy = 0;
    
    for(size_t i = 0 ;  i < N ; ++i)
        kineticEnergy += 0.5*mass[i]*vel[i].normSquare();
    
    return kineticEnergy;
}

template<typename Scalar, size_t N>
double getPotentialEnergy(const std::array<Scalar,N>& mass, const std::array<vec3<Scalar>, N>& pos)
{
    Scalar potentialEnergy = 0;
    
    for(size_t i = 0 ;  i < N ; ++i)
        for(size_t j = i + 1; j < N; ++j)
            potentialEnergy -= mass[i]*mass[j]/distance(pos[i], pos[j]);
    
    return potentialEnergy;
}
template<typename Scalar, size_t N>
inline double getTotalEnergy(const std::array<Scalar,N>& mass, const std::array<vec3<Scalar>, N>& pos, const std::array<vec3<Scalar>, N>& vel)
{
    Scalar potentialEnergy = 0;
    Scalar kineticEnergy   = 0;
    
    for(size_t i = 0 ;  i < N ; ++i)
    {
        kineticEnergy += 0.5*mass[i]*vel[i].normSquare();
        
        for(size_t j = i + 1; j < N; ++j)
            potentialEnergy -= mass[i]*mass[j]/distance(pos[i], pos[j]);
    }
    return potentialEnergy + kineticEnergy;
}

template<typename T>
void print(T& var)
{
    for(size_t i = 0 ; i < var.size(); ++i )
        std::cout << var[i] << ' ';
    std::cout << '\n';
}
#endif
