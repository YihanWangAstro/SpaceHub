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

/**
 *  @brief Class of dynamical variable.
 *
 *  All variables in this class are physical quantities need to be evolved. If you want to create you own dynamic state,
 *  make sure remove constant quantities away from the class. The interface array() can be used to operate the data in
 *  the class as a plain array(for evolution). Extra constant physical quantities waste the calculation.
 */
template <typename DataType, size_t N>
class dynamics
{
public:
    typedef DataType                   Scalar;
    typedef vec3<Scalar>               Vector;
    typedef std::array<vec3<Scalar>,N> VectorArray;
    typedef std::array<Scalar,N>       ScalarArray;
    
    /** @brief Get the number of the particles.
     *  @return The particle number.
     */
    constexpr static size_t size(){return N;}
    
    /** @brief Get the total data number.
     *  @return The data number.
     */
    constexpr static size_t volume(){return 6*N + 1;}
    
    /** @brief Array of position of the particles. Element is 3D vector.*/
    VectorArray pos;
    
    /** @brief Array of velocity of the particles. Element is 3D vector.*/
    VectorArray vel;
    
    /** @brief The physical time of the dynamic system*/
    Scalar      time{0.0};

    /** @brief Transfer this class to a plain array.
     *  @return The reference of head of this class, reinterpret as a plain array.
     */
    std::array<Scalar, volume()>& array()
    {
        return reinterpret_cast<std::array<Scalar, volume()>&>(*this);
    }
    
    /** @brief Initialize extra user defined variables. Interface required for other class.
     *  @param mass The mass of particles, might be required for initialization.
     */
    void initAddiVariable(ScalarArray& mass){};
    
    /** @brief Set all data to be zero.*/
    void setZero()
    {
        memset(&pos[0], 0, sizeof(Scalar)*volume());
    }
};
#endif

