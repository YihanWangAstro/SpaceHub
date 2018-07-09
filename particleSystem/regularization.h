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
    using Scalar      = DynamicState::Scalar;
    using Vector      = DynamicState::Vector;
    using VectorArray = DynamicState::VectorArray;
    using ScalarArray = DynamicState::ScalarArray;
    using IntArray    = DynamicState::IntArray;

    /** @brief Calculate the physical time for position advance from integration step size
     *  @param mass     Array of particle mass.
     *  @param dyn      Dynamic system contains position, velocity and regularization variables.
     *                  See example class in dynamicState.h.
     *  @param stepSize Integration step size. This could not be the physical time. Look
     *                  references for details in class despriction.
     */
    Scalar getPhysicalPosTime(DynamicState& dyn, Scalar stepSize)
    {
        return stepSize / (dyn.bindE() + getKineticEnergy(dyn.mass(), dyn.vel()));
    }

    /** @brief Calculate the physical time for velocity advance from integration step size
     *  @param mass     Array of particle mass.
     *  @param dyn      Dynamic system contains position, velocity and regularization variables.
     *                  See example class in dynamicState.h.
     *  @param stepSize Integration step size. This could not be the physical time. Look
     *                  references for details in class despriction.
     */
    Scalar getPhysicalVelTime(DynamicState& dyn, Scalar stepSize)
    {
        return stepSize / (-getPotentialEnergy(dyn.mass(), dyn.pos()));
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
    using Scalar      = DynamicState::Scalar;
    using Vector      = DynamicState::Vector;
    using VectorArray = DynamicState::VectorArray;
    using ScalarArray = DynamicState::ScalarArray;
    using IntArray    = DynamicState::IntArray;

    /** @brief Calculate the physical time for position advance from integration step size
     *  @param mass     Array of particle mass.
     *  @param dyn      Dynamic system contains position, velocity and regularization variables.
     *                  See example class in dynamicState.h.
     *  @param stepSize Integration step size. This could not be the physical time. Look
     *                  references for details in class despriction.
     */
    Scalar getPhysicalPosTime(DynamicState& dyn, Scalar stepSize)
    {
        return stepSize / dyn.omega();
    }

    /** @brief Calculate the physical time for velocity advance from integration step size
     *  @param mass     Array of particle mass.
     *  @param dyn      Dynamic system contains position, velocity and regularization variables.
     *                  See example class in dynamicState.h.
     *  @param stepSize Integration step size. This could not be the physical time. Look
     *                  references for details in class despriction.
     */
    Scalar getPhysicalVelTime(DynamicState& dyn, Scalar stepSize)
    {
        return stepSize / dyn.capitalOmega();
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
    using Scalar      = DynamicState::Scalar;
    using Vector      = DynamicState::Vector;
    using VectorArray = DynamicState::VectorArray;
    using ScalarArray = DynamicState::ScalarArray;
    using IntArray    = DynamicState::IntArray;

    /** @brief Calculate the physical time for position advance from integration step size
     *  @param mass     Array of particle mass.
     *  @param dyn      Dynamic system contains position, velocity and regularization variables.
     *                  See example class in dynamicState.h.
     *  @param stepSize Integration step size.
     */
    Scalar getPhysicalPosTime(DynamicState& dyn, Scalar stepSize)
    {
        return stepSize;
    }

    /** @brief Calculate the physical time for velocity advance from integration step size
     *  @param mass     Array of particle mass.
     *  @param dyn      Dynamic system contains position, velocity and regularization variables.
     *                  See example class in dynamicState.h.
     *  @param stepSize Integration step size.
     */
    Scalar getPhysicalVelTime(DynamicState& dyn, Scalar stepSize)
    {
        return stepSize;
    }
};
#endif

