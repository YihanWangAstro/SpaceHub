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
    typedef DataType                   Scalar;
    typedef vec3<Scalar>               Vector;
    typedef std::array<vec3<Scalar>,N> VectorArray;
    typedef std::array<Scalar,N>       ScalarArray;
    typedef std::array<size_t,N>       IndexArray;
    
    /** @brief Get the number of the particles.
     *  @return The particle number.
     */
    constexpr static size_t size(){return N;}
    
    /** @brief Get the total data number.
     *  @return The data number.
     */
    constexpr static size_t volume(){return 6*N + 3;}
    
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
     *  @return The reference of head of this class, reinterpret as a plain array.
     */
    std::array<Scalar, volume()>& array()
    {
        return reinterpret_cast<std::array<Scalar, volume()>&>(pos);
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
        bindE = -getTotalEnergy(mass,pos,vel);
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

