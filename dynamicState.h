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
    typedef DataType                      Scalar;
    typedef vec3<Scalar>                  Vector;
    typedef std::array<vec3<Scalar>, N>   VectorArray;
    typedef std::array<Scalar, N>         ScalarArray;
    typedef std::array<Scalar, 6 * N + 1> PlainArray;
    /** @brief Get the number of the particles.
     *  @return The particle number.
     */
    constexpr static size_t size()
    {
        return N;
    }

    /** @brief Get the total data number.
     *  @return The data number.
     */
    constexpr static size_t volume()
    {
        return 6 * N + 1;
    }

    /** @brief Array of position of the particles. Element is 3D vector.*/
    VectorArray pos;

    /** @brief Array of velocity of the particles. Element is 3D vector.*/
    VectorArray vel;

    /** @brief The physical time of the dynamic system*/
    Scalar      time{0.0};

    /** @brief Transfer this class to a plain array.
     *  @param arr The destination plain array.
     */
    void flatten(PlainArray& arr)
    {
        size_t len = sizeof(VectorArray);
        memcpy(static_cast<void*>(&arr[0]), static_cast<void*>(&pos[0]), len);
        memcpy(static_cast<void*>(&arr[3*N]), static_cast<void*>(&vel[0]), len);
        memcpy(static_cast<void*>(&arr[6*N]), static_cast<void*>(&time), sizeof(Scalar));
    }
    
    /** @brief Load data from a plain array.
     *  @param arr The plain array data.
     */
    void loadFlatten(PlainArray& arr)
    {
        size_t len = sizeof(VectorArray);
        memcpy(static_cast<void*>(&pos[0]), static_cast<void*>(&arr[0]), len);
        memcpy(static_cast<void*>(&vel[0]), static_cast<void*>(&arr[3*N]),  len);
        memcpy(static_cast<void*>(&time), static_cast<void*>(&arr[6*N]), sizeof(Scalar));
    }
    

    /** @brief Initialize extra user defined variables. Interface required for other class.
     *  @param mass The mass of particles, might be required for initialization.
     */
    void initAddiVariable(ScalarArray& mass) {};

    /** @brief Set all data to be zero.*/
    void setZero()
    {
        memset(&pos[0], 0, sizeof(Scalar)*volume());
    }
};
#endif

