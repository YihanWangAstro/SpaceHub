////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Filename:vector.h                                                                                                   //
//Author:Yihan Wang                                                                                                   //
//                                                                                                                    //
//                                                                                                                    //
//Description:                                                                                                        //
//  Basic 3-d vector class for all vector operations.                                                                 //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef REGULARSTATE_H
#define REGULARSTATE_H
#include "chain.h"
#include "../libs.h"

/**
 *  @brief Class of dynamical system with regularization variables.
 *
 *  A simple extension of class dynamics in dynamicState.h. Used for regularization system. See detail
 *  in https://academic.oup.com/mnras/article/372/1/219/974304 .
 */
template <typename DataType, size_t N>
class reguDynamics
{
public:
    typedef DataType                      Scalar;
    typedef vec3<Scalar>                  Vector;
    typedef std::array<vec3<Scalar>, N>   VectorArray;
    typedef std::array<Scalar, N>         ScalarArray;
    typedef std::array<size_t, N>         IndexArray;
    typedef std::array<Scalar, 6 * N + 3> PlainArray;
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
        return 6 * N + 3;
    }

    /** @brief Array of position of the particles. Element is 3D vector.*/
    VectorArray pos;

    /** @brief Array of velocity of the particles. Element is 3D vector.*/
    VectorArray vel;

    /** @brief The physical time of the dynamic system*/
    Scalar time{0.0};

    /** @brief The binding energy(for regularization) of the dynamic system*/
    Scalar bindE{0.0};

    /** @brief The regularization variable of the dynamic system*/
    Scalar omega{0.0};

    /** @brief Transfer this class to a plain array.
     *  @param arr The destination plain array.
     */
    void flatten(PlainArray& arr)
    {
        size_t len = sizeof(VectorArray);
        memcpy(static_cast<void*>(&arr[0]), static_cast<void*>(&pos[0]), len);
        memcpy(static_cast<void*>(&arr[3*N]), static_cast<void*>(&vel[0]), len);
        memcpy(static_cast<void*>(&arr[6*N]), static_cast<void*>(&time), sizeof(Scalar));
        memcpy(static_cast<void*>(&arr[6*N + 1]), static_cast<void*>(&bindE), sizeof(Scalar));
        memcpy(static_cast<void*>(&arr[6*N + 2]), static_cast<void*>(&omega), sizeof(Scalar));
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
        memcpy(static_cast<void*>(&bindE), static_cast<void*>(&arr[6*N + 1]), sizeof(Scalar));
        memcpy(static_cast<void*>(&omega), static_cast<void*>(&arr[6*N + 2 ]), sizeof(Scalar));
    }

    /** @brief Set all data to be zero.*/
    void setZero()
    {
        memset(&pos[0], 0, sizeof(Scalar)*volume());
    }

    /** @brief Calculate the regularization variable omega
     *
     *  @param mass The mass of particles, might be required for calculation.
     *  @return The calculated value of omega.
     */
    Scalar getOmega(ScalarArray& mass)
    {
        return -getPotentialEnergy(mass, pos);
    }

    /** @brief Initialize extra user defined variables. Interface required for other class.
     *
     *  Initialize regularizaiton variable bindE and omega.
     *  @param mass The mass of particles, might be required for initialization.
     */
    void initAddiVariable(ScalarArray& mass)
    {
        bindE = -getTotalEnergy(mass, pos, vel);
        omega = getOmega(mass);
    }

    /** @brief Transfer Cartesian coordinate regularization system to chain regularization system.
     *
     *  Coordinate transformation. From Cartesian to chain. See details in
     *  https://link.springer.com/article/10.1007%2FBF00695714 .
     *  @param chainData The destination regularization system in chain coordinates.
     *  @param index     The maping index between Cartesian coordinates and chain coordinates.
     */
    void toChain(reguDynamics& chainData, IndexArray& index)
    {
        chain::synChain(pos, chainData.pos, index);
        chain::synChain(vel, chainData.vel, index);
        chainData.time  = time;
        chainData.omega = omega;
        chainData.bindE = bindE;
    }

    /** @brief Transfer chain coordinate regularization system to Cartesian regularization system.
     *
     *  Coordinate transformation. From chain to Cartesian. See details in
     *  https://link.springer.com/article/10.1007%2FBF00695714 .
     *  @param cartesian The destination regularization system in Cartesian coordinates.
     *  @param index     The maping index between Cartesian coordinates and chain coordinates.
     */
    void toCartesian(reguDynamics& cartesian, IndexArray& index)
    {
        chain::synCartesian(pos, cartesian.pos, index);
        chain::synCartesian(vel, cartesian.vel, index);
        cartesian.time  = time;
        cartesian.omega = omega;
        cartesian.bindE = bindE;
    }

    /** @brief Move particles to central mass coordinates
     *
     *  Move position and velocity to central mass coordinates.
     *  @param mass Mass of the particles required for moving.
     */
    void moveToCentralMassCoords(ScalarArray& mass)
    {
        MoveToCentralMassCoordinate(mass, pos);
        MoveToCentralMassCoordinate(mass, vel);
    }
};
#endif

