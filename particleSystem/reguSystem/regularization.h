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
#include "../../libs.h"

namespace SpaceH
{
/**  @brief logH extention algorithmatic regularization interface
  *
  *  See detials in https://link.springer.com/article/10.1023%2FA%3A1008368322547 and
  *                 http://iopscience.iop.org/article/10.1086/301102/meta .
 */
template<typename Particles>
class LogH
{
public:
    /* Typedef */
    using type   = typename Particles::type;
    using Scalar = typename type::Scalar;
    /* Typedef */

    /** @brief Calculate the physical time for position advance from integration step size
     *  @param partc    Dynamic system contains position, velocity and regularization variables.
     *                  See example class in particles.h.
     *  @param stepSize Integration step size. This could not be the physical time. Look
     *                  references for details in class despriction.
     */
    inline Scalar getPhysicalPosTime(const Particles& partc, Scalar stepSize)
    {
        return stepSize / (partc.bindE() + getKineticEnergy(partc.mass(), partc.vel()));
    }

    /** @brief Calculate the physical time for velocity advance from integration step size
     *  @param partc    Dynamic system contains position, velocity and regularization variables.
     *                  See example class in particles.h.
     *  @param stepSize Integration step size. This could not be the physical time. Look
     *                  references for details in class despriction.
     */
    inline Scalar getPhysicalVelTime(const Particles& partc, Scalar stepSize)
    {
        return stepSize / (-getPotentialEnergy(partc.mass(), partc.pos()));
    }
};

/**  @brief Time Transform Leapfrog algorithmatic regularization interface
 *
 *  See detials in https://link.springer.com/article/10.1023%2FA%3A1021149313347 .
 */
template<typename Particles>
class TTL
{
public:
    /* Typedef */
    using type   = typename Particles::type;
    using Scalar = typename type::Scalar;
    /* Typedef */

    /** @brief Calculate the physical time for position advance from integration step size
     *  @param partc    Dynamic system contains position, velocity and regularization variables.
     *                  See example class in particles.h.
     *  @param stepSize Integration step size. This could not be the physical time. Look
     *                  references for details in class despriction.
     */
    inline Scalar getPhysicalPosTime(const Particles& partc, Scalar stepSize)
    {
        return stepSize / partc.omega();
    }

    /** @brief Calculate the physical time for velocity advance from integration step size
     *  @param partc    Dynamic system contains position, velocity and regularization variables.
     *                  See example class in particles.h.
     *  @param stepSize Integration step size. This could not be the physical time. Look
     *                  references for details in class despriction.
     */
    inline Scalar getPhysicalVelTime(const Particles& partc, Scalar stepSize)
    {
        return stepSize / partc.getCapitalOmega();
    }
};

/**  @brief Ordinary algorithmatic regularization interface
 *
 *   No regularization.
 */
template<typename Particles>
class NoRegu
{
public:
    /* Typedef */
    using type   = typename Particles::type;
    using Scalar = typename type::Scalar;
    /* Typedef */

    /** @brief Calculate the physical time for position advance from integration step size
     *  @param partc    Dynamic system contains position, velocity and regularization variables.
     *                  See example class in particles.h.
     *  @param stepSize Integration step size.
     */
    inline Scalar getPhysicalPosTime(const Particles& partc, Scalar stepSize)
    {
        return stepSize;
    }

    /** @brief Calculate the physical time for velocity advance from integration step size
     *  @param partc    Dynamic system contains position, velocity and regularization variables.
     *                  See example class in particles.h.
     *  @param stepSize Integration step size.
     */
    inline Scalar getPhysicalVelTime(const Particles& partc, Scalar stepSize)
    {
        return stepSize;
    }
};
}
#endif

