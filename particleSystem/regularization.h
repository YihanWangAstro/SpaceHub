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

/**  @brief logH extention algorithmatic regularization interface
  *
  *  See detials in https://link.springer.com/article/10.1023%2FA%3A1008368322547 and
  *                 http://iopscience.iop.org/article/10.1086/301102/meta .
 */
template<typename DynamicState>
class logH
{
public:
    typedef typename DynamicState::Scalar Scalar;
    constexpr static size_t size(){return DynamicState::size();}
    
    /** @brief Calculate the physical time for position advance from integration step size
     *  @param mass     Array of particle mass.
     *  @param dyn      Dynamic system contains position, velocity and regularization variables.
     *                  See example class in dynamicState.h.
     *  @param stepSize Integration step size. This could not be the physical time. Look
     *                  references for details in class despriction.
     */
    Scalar getPhysicalPosTime(std::array<Scalar, size()>& mass, DynamicState& dyn, Scalar stepSize)
    {
        return stepSize/(dyn.bindE + getKineticEnergy(mass, dyn.vel));
    }
    
    /** @brief Calculate the physical time for velocity advance from integration step size
     *  @param mass     Array of particle mass.
     *  @param dyn      Dynamic system contains position, velocity and regularization variables.
     *                  See example class in dynamicState.h.
     *  @param stepSize Integration step size. This could not be the physical time. Look
     *                  references for details in class despriction.
     */
    Scalar getPhysicalVelTime(std::array<Scalar, size()>& mass, DynamicState& dyn, Scalar stepSize)
    {
        return stepSize/(-getPotentialEnergy(mass, dyn.pos));
    }
};

/**  @brief Time Transform Leapfrog algorithmatic regularization interface
 *
 *  See detials in https://link.springer.com/article/10.1023%2FA%3A1021149313347 .
 */
template<typename DynamicState>
class TTL
{
public:
    typedef typename DynamicState::Scalar Scalar;
    constexpr static size_t size(){return DynamicState::size();}
    
    /** @brief Calculate the physical time for position advance from integration step size
     *  @param mass     Array of particle mass.
     *  @param dyn      Dynamic system contains position, velocity and regularization variables.
     *                  See example class in dynamicState.h.
     *  @param stepSize Integration step size. This could not be the physical time. Look
     *                  references for details in class despriction.
     */
    Scalar getPhysicalPosTime(std::array<Scalar, size()>& mass, DynamicState& dyn, Scalar stepSize)
    {
        return stepSize/dyn.omega;
    }
    
    /** @brief Calculate the physical time for velocity advance from integration step size
     *  @param mass     Array of particle mass.
     *  @param dyn      Dynamic system contains position, velocity and regularization variables.
     *                  See example class in dynamicState.h.
     *  @param stepSize Integration step size. This could not be the physical time. Look
     *                  references for details in class despriction.
     */
    Scalar getPhysicalVelTime(std::array<Scalar, size()>& mass, DynamicState& dyn, Scalar stepSize)
    {
        return stepSize/dyn.getOmega(mass);
    }
};

/**  @brief Ordinary algorithmatic regularization interface
 *
 *   No regularization.
 */
template<typename DynamicState>
class NoRegu
{
public:
    typedef typename DynamicState::Scalar Scalar;
    constexpr static size_t size(){return DynamicState::size();}
    
    /** @brief Calculate the physical time for position advance from integration step size
     *  @param mass     Array of particle mass.
     *  @param dyn      Dynamic system contains position, velocity and regularization variables.
     *                  See example class in dynamicState.h.
     *  @param stepSize Integration step size.
     */
    Scalar getPhysicalPosTime(std::array<Scalar, size()>& mass, DynamicState& dyn, Scalar stepSize)
    {
        return stepSize;
    }
    
    /** @brief Calculate the physical time for velocity advance from integration step size
     *  @param mass     Array of particle mass.
     *  @param dyn      Dynamic system contains position, velocity and regularization variables.
     *                  See example class in dynamicState.h.
     *  @param stepSize Integration step size. 
     */
    Scalar getPhysicalVelTime(std::array<Scalar, size()>& mass, DynamicState& dyn, Scalar stepSize)
    {
        return stepSize;
    }
};
#endif

