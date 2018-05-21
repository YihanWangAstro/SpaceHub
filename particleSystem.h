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

/**  @brief Base class of particle System.
 *
 *   Base particles system class. Other particle system can inherit this class. Considering the performance, we don't
 *   set virtual function. Here we use CRTP technique to bind derived class method. See more details in
 *   https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern .
 */
template <typename Derived, typename EvolvedData>
class particleSystem
{
public:

    typedef typename EvolvedData::Scalar              Scalar;
    typedef typename EvolvedData::Vector              Vector;
    typedef typename EvolvedData::VectorArray         VectorArray;
    typedef typename EvolvedData::ScalarArray         ScalarArray;
    typedef std::array<size_t, EvolvedData::size()>   IntArray;
    typedef std::array<Scalar, EvolvedData::volume()> PlainArray;

    /**  @brief Evolved variables class*/
    EvolvedData  dynState;

    /**  @brief Position array interface. Reference to dynState.pos*/
    VectorArray& pos;

    /**  @brief Velocity array interface. Reference to dynState.vel*/
    VectorArray& vel;

    /**  @brief Physical time scalar interface. Reference to dynState.time*/
    Scalar&      time;

    /**  @brief Mass array.*/
    ScalarArray  mass;

    /**  @brief Radius array.*/
    ScalarArray  radius;

    /**  @brief Particle type array.*/
    IntArray     type;

    /**  @brief Default construction.
     *
     *   Default constructor. Bind position, velocity and time interface reference to dynState.pos, dynState.vel and
     *   dyn.time.
     */
    particleSystem() : pos(dynState.pos), vel(dynState.vel), time(dynState.time) {}

    /** @brief Virtualize default destructor.*/
    virtual ~particleSystem() {}

    /** @brief Get the number of the particles.
     *  @return The particle number.
     */
    constexpr static size_t size()
    {
        return EvolvedData::size();
    }

    /** @brief Get the total dynamic scalar number.
     *  @return The dynamic scalar number.
     */
    constexpr static size_t volume()
    {
        return EvolvedData::volume();
    }

    /** @brief Transfer evolved data to a plain array.
     *  @return The reference of a plain array.
     */
    PlainArray&   array()
    {
        return dynState.array();
    }

    /** @brief Interface to rescale the time.
     *
     *  Interace used by dynamic system. Transfer integration time(For some system, integration time is different from
     *  physical time) to physical time.
     *  @return The phsyical time.
     */
    Scalar timeScale(Scalar scale)
    {
        return scale;
    }

    /** @brief Load data from a plain array.*/
    void load(PlainArray& data);

    /** @brief Output data to standard c++ ostream. */
    std::ostream& write(std::ostream&)const;

    /** @brief Input data from standard c++ istream.*/
    std::istream& read (std::istream&);

    /** @brief Overload operator = */
    const particleSystem& operator=(const particleSystem& other);
protected:
    /**  @brief Acceleration array used to update velocity.*/
    VectorArray  acc;
};

/** @brief Overload operator << */
template <typename Derived, typename EvolvedData>
std::ostream& operator<<(std::ostream& os, const particleSystem<Derived, EvolvedData>& data)
{
    return static_cast<const Derived*>(&data)->write(os);
}

/** @brief Overload operator >> */
template <typename Derived, typename EvolvedData>
std::istream& operator>>(std::istream& is, particleSystem<Derived, EvolvedData>& data)
{
    return static_cast<Derived*>(&data)->read(is);
}

/** @brief Overload operator = */
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

/** @brief Input data from standard c++ istream.
 *
 *  Implement of CRTP '>>' method.
 */
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

/** @brief Output data to standard c++ ostream.
 *
 *  Implement of CRTP '<<' method.
 */
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

/** @brief Load data from a plain array.
 *
 *  Interface usded by integrator and ODE iterator. Load data from a plain array processed by itegrator and
 *  iterator. Derived class could overload this function to additional process.
 *
 *  @param data Plain scalar array.
 */
template <typename Derived, typename EvolvedData>
void particleSystem<Derived, EvolvedData>::load(PlainArray& data)
{
    dynState.array() = data;
}
#endif
