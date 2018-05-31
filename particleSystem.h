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
 *   set virtual function.
 */
template <typename EvolvedData>
class particleSystem
{
public:

    typedef typename EvolvedData::Scalar            Scalar;
    typedef typename EvolvedData::Vector            Vector;
    typedef typename EvolvedData::VectorArray       VectorArray;
    typedef typename EvolvedData::ScalarArray       ScalarArray;
    typedef typename EvolvedData::PlainArray        PlainArray;
    typedef std::array<size_t, EvolvedData::size()> IndexArray;

    /**  @brief Physical time scalar interface. Reference to dynState.time*/
    inline Scalar& time(){ return dynState.time; }
    
    /**  @brief Position array interface. Reference to dynState.pos*/
    inline VectorArray& pos(){ return dynState.pos; }
    
    /**  @brief Velocity array interface. Reference to dynState.vel*/
    inline VectorArray& vel(){ return dynState.vel; }
    
    /**  @brief Mass array interface. Reference to m.*/
    inline ScalarArray& mass(){ return m; }
    
    /**  @brief Radius array interface. Reference to r.*/
    inline ScalarArray& radius(){ return rad; }
    
    /**  @brief Particle type array interface. Reference to t.*/
    inline IndexArray& type(){ return tp; }
    
    /**  @brief Position vector interface. Reference to dynState.pos[i]*/
    inline Vector& pos(size_t i){ return dynState.pos[i]; }
    
    /**  @brief Velocity vecotr interface. Reference to dynState.vel[i]*/
    inline Vector& vel(size_t i){ return dynState.vel[i]; }
    
    /**  @brief Mass interface. Reference to m.*/
    inline Scalar& mass(size_t i){ return m[i]; }
    
    /**  @brief Radius interface. Reference to r.*/
    inline Scalar& radius(size_t i){ return rad[i]; }
    
    /**  @brief Particle type interface. Reference to t.*/
    inline size_t& type(size_t i){ return tp[i]; }
    
    /**  @brief Physical time scalar const interface. Reference to dynState.time*/
    inline const Scalar& time() const { return dynState.time; }
    
    /**  @brief Position array const interface. Reference to dynState.pos*/
    inline const VectorArray& pos() const { return dynState.pos; }
    
    /**  @brief Velocity array const interface. Reference to dynState.vel*/
    inline const VectorArray& vel() const { return dynState.vel; }
    
    /**  @brief Mass array const interface. Reference to m.*/
    inline const ScalarArray& mass() const { return m; }
    
    /**  @brief Radius array const interface. Reference to r.*/
    inline const ScalarArray& radius() const { return rad; }
    
    /**  @brief Particle type array const interface. Reference to t.*/
    inline const IndexArray& type() const { return tp; }
    
    /**  @brief Position vector const interface. Reference to dynState.pos[i]*/
    inline const Vector& pos(size_t i) const { return dynState.pos[i]; }
    
    /**  @brief Velocity vecotr const interface. Reference to dynState.vel[i]*/
    inline const Vector& vel(size_t i) const { return dynState.vel[i]; }
    
    /**  @brief Mass const interface. Reference to m.*/
    inline const  Scalar& mass(size_t i) const { return m[i]; }
    
    /**  @brief Radius const interface. Reference to r.*/
    inline const Scalar& radius(size_t i) const { return rad[i]; }
    
    /**  @brief Particle type const interface. Reference to t.*/
    inline const size_t& type(size_t i) const { return tp[i]; }

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
     *  @param arr The destination plain array.
     */
    void flatten(PlainArray& arr)
    {
        dynState.flatten(arr);
    }
    
    /** @brief Load data from a plain array.
     *  @param arr The plain array data.
     */
    void loadFlatten(PlainArray& arr)
    {
         dynState.loadFlatten(arr);
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
    
    /** @brief Preprocess before iteration*/
    void preIterProcess(){}
    
    /** @brief After process after iteration*/
    void afterIterProcess(){}
    
    /** @brief Output data to standard c++ ostream. */
    std::ostream& write(std::ostream&)const;
    
    /** @brief Input data from standard c++ istream.*/
    std::istream& read (std::istream&);
    
    /** @brief Virtualize default destructor.*/
    virtual ~particleSystem() {}
    
    /** @brief Overload operator << */
    friend std::ostream& operator<<(std::ostream& output, const particleSystem& sys)
    {
        return sys.write(output);
    }
    /** @brief Input from istream */
    friend std::istream& operator>>(std::istream& input, particleSystem& sys)
    {
        return sys.read(input);
    }
    
protected:
    /**  @brief Evolved variables class*/
    EvolvedData dynState;
    
    /**  @brief Acceleration array used to update velocity.*/
    VectorArray acc;
    
    /**  @brief Mass array.*/
    ScalarArray m;
    
    /**  @brief Radius array.*/
    ScalarArray rad;
    
    /**  @brief Particle type array.*/
    IndexArray tp;
};

/** @brief Input data from standard c++ istream.
 *
 *  Implement of CRTP '>>' method.
 */
template <typename EvolvedData>
std::istream& particleSystem<EvolvedData>::read(std::istream& input)
{
    input >> dynState.time;
    size_t id;

    for(size_t i = 0 ; i < size() ; ++i)
        input >> id >> tp[i] >> m[i] >> rad[i] >> dynState.pos[i] >> dynState.vel[i];

    memset(&(this->acc[0]), 0, sizeof(Vector)*size());
    dynState.moveToCentralMassCoords(m);
    dynState.initAddiVariable(m);
    return input;
}

/** @brief Output data to standard c++ ostream.
 *
 *  Implement of CRTP '<<' method.
 */
template <typename EvolvedData>
std::ostream& particleSystem<EvolvedData>::write(std::ostream& output) const
{
    output << "#" << size() << " " << time() << "\r\n";

    for(size_t i = 0 ; i < size() ; ++i)
    {
        output << i << " " << tp[i] << " " << m[i] << " " << rad[i]
                    << " " << dynState.pos[i]  << " " << dynState.vel[i]  << "\r\n";
    }

    return output;
}
#endif
