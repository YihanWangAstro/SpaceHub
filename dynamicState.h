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

template <typename DataType, size_t N>
class dynamics
{
public:
    typedef DataType                   Scalar;
    typedef vec3<Scalar>               Vector;
    typedef std::array<vec3<Scalar>,N> VectorArray;
    typedef std::array<Scalar,N>       ScalarArray;
    
    constexpr static size_t size(){return N;}
    constexpr static size_t volume(){return 6*N + 1;}
    
    VectorArray pos;
    VectorArray vel;
    Scalar      time{0.0};
    
    dynamics(){}
    dynamics(const dynamics& ds) : pos(ds.pos), vel(ds.vel), time(ds.time){};
    std::array<Scalar, volume()>& array(){ return reinterpret_cast<std::array<Scalar, volume()>&>(*this);}
    void initAddiVariable(ScalarArray& mass){};
    void setZero(){ memset(&pos[0], 0, sizeof(Scalar)*volume());}
};

template <typename DataType, size_t N>
class reguDynamics
{
public:
    typedef DataType                   Scalar;
    typedef vec3<Scalar>               Vector;
    typedef std::array<vec3<Scalar>,N> VectorArray;
    typedef std::array<Scalar,N>       ScalarArray;
    
    constexpr static size_t size(){return N;}
    constexpr static size_t volume(){return 6*N + 3;}
    
    VectorArray pos;
    VectorArray vel;
    Scalar      time{0.0};
    Scalar      bindE{0.0};
    Scalar      omega{0.0};
    
    reguDynamics(){}
    reguDynamics(const reguDynamics& ds) : pos(ds.pos), vel(ds.vel), time(ds.time), bindE(ds.bindE), omega(ds.omega){};
    std::array<Scalar, volume()>& array(){ return reinterpret_cast<std::array<Scalar, volume()>&>(pos);}
    void setZero(){memset(&pos[0], 0, sizeof(Scalar)*volume());}
    Scalar getOmega(ScalarArray& mass)
    {
        return -getPotentialEnergy(mass, pos);
    }
    void initAddiVariable(ScalarArray& mass)
    {
        bindE = -getTotalEnergy(mass,pos,vel);
        omega = getOmega(mass);
    }
};

template <typename DataType, size_t N>
class GAR
{
public:
    typedef DataType                   Scalar;
    typedef vec3<Scalar>               Vector;
    typedef std::array<vec3<Scalar>,N> VectorArray;
    typedef std::array<Scalar,N>       ScalarArray;
    
    constexpr static size_t size(){return N;}
    constexpr static size_t volume(){return 9*N + 3;}
    
    VectorArray pos;
    VectorArray vel;
    VectorArray auxiVel;
    Scalar      time{0.0};
    Scalar      bindE{0.0};
    Scalar      omega{0.0};
    
    GAR(){}
    GAR(const GAR& ds) : pos(ds.pos), vel(ds.vel), time(ds.time), bindE(ds.bindE), omega(ds.omega){};
    std::array<Scalar, volume()>& array(){ return reinterpret_cast<std::array<Scalar, volume()>&>(pos);}
    void setZero(){memset(&pos[0], 0, sizeof(Scalar)*volume());}
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
};
#endif

