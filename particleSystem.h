////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Filename:particleSystem.h                                                                                           //
//Author:Yihan Wang                                                                                                   //
//                                                                                                                    //
//                                                                                                                    //
//Description:                                                                                                        //
//  The main class of this n-body code. This class includes                                                           //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef GENPARTICLESYSTEM_H
#define GENPARTICLESYSTEM_H
#include "vector3.h"
#include "libs.h"
#include "errhand.h"
#include <array>
#include <fstream>
#include <cstring>
#include <iomanip>

template <typename Derived, typename EvolvedData>
class particleSystem
{
public:
    /////////////////////////////////Type Define////////////////////////////////////
    typedef typename EvolvedData::Scalar              Scalar;
    typedef typename EvolvedData::Vector              Vector;
    typedef typename EvolvedData::VectorArray         VectorArray;
    typedef typename EvolvedData::ScalarArray         ScalarArray;
    typedef std::array<size_t, EvolvedData::size()>    IntArray;
    typedef std::array<Scalar, EvolvedData::volume()> PlainArray;
    ///////////////////////////////////Interface////////////////////////////////////
    constexpr static size_t size()
    {
        return EvolvedData::size();
    }
    constexpr static size_t volume()
    {
        return EvolvedData::volume();
    }
    EvolvedData  dynState;
    VectorArray& pos;
    VectorArray& vel;
    Scalar&      time;
protected:
    VectorArray  acc;
public:
    ScalarArray  mass;
    ScalarArray  radius;
    IntArray     type;
    particleSystem() : pos(dynState.pos), vel(dynState.vel), time(dynState.time) {}
    virtual ~particleSystem() {}
    
    std::ostream& write(std::ostream&)const;
    std::istream& read (std::istream&);
    PlainArray&   array()
    {
        return dynState.array();
    }
    Scalar timeScale(Scalar scale);
    void load(PlainArray& data);
    const particleSystem& operator=(const particleSystem& other);
};

//////////////////////////////CRTP Interface////////////////////////////////////
template <typename Derived, typename EvolvedData>
std::ostream& operator<<(std::ostream& os, const particleSystem<Derived, EvolvedData>& data)
{
    return static_cast<const Derived*>(&data)->write(os);
}

template <typename Derived, typename EvolvedData>
std::istream& operator>>(std::istream& is, particleSystem<Derived, EvolvedData>& data)
{
    return static_cast<Derived*>(&data)->read(is);
}

//////////////////////////////Default implement/////////////////////////////////
template <typename Derived, typename EvolvedData>
const particleSystem<Derived, EvolvedData>& particleSystem<Derived, EvolvedData>::operator=(const particleSystem& other)
{
    dynState = other.dynState;
    acc      = other.acc;
    mass     = other.mass;
    radius   = other.radius;
    type     = other.type;
    return *this;
}

template <typename Derived, typename EvolvedData>
typename EvolvedData::Scalar particleSystem<Derived, EvolvedData>::timeScale(Scalar scale)
{
    return scale;
}

template <typename Derived, typename EvolvedData>
std::istream& particleSystem<Derived, EvolvedData>::read(std::istream& input)
{
    input >> time;
    size_t id;
    
    for(size_t i = 0 ; i < size() ; ++i)
        input >> id >> type[i] >> mass[i] >> radius[i] >> pos[i] >> vel[i];
        
    memset(acc, 0, sizeof(Scalar) * 3 * size());
    MoveToCentralMassCoordinate(mass, dynState.pos);
    MoveToCentralMassCoordinate(mass, dynState.vel);
    dynState.initAddiVariable  (mass);
    return input;
}

template <typename Derived, typename EvolvedData>
std::ostream& particleSystem<Derived, EvolvedData>::write(std::ostream& output) const
{
    output << "#" << size() << " " << time << "\r\n";
    
    for(size_t i = 0 ; i < size() ; ++i)
    {
        output << i << " " << type[i] << " " << mass[i] << " " << radius[i]
               << " " << pos [i] << " " << vel [i] << "\r\n";
    }
    
    return output;
}

template <typename Derived, typename EvolvedData>
void particleSystem<Derived, EvolvedData>::load(PlainArray& data)
{
    dynState.array() = data;
}
#endif
