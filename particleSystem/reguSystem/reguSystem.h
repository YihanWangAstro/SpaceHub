////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Filename:chainSystem.h                                                                                              //
//Author:Yihan Wang                                                                                                   //
//                                                                                                                    //
//                                                                                                                    //
//Description:                                                                                                        //
//  The main class of this n-body code. This class includes                                                           //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef REGULARIZEDSYSTEM_H
#define REGULARIZEDSYSTEM_H
#include "../../particleSystem.h"


/**  @brief Regularized particle System.
 *
 *   Regularied particle system. See details in https://link.springer.com/article/10.1023%2FA%3A1008368322547 ,
 *   http://iopscience.iop.org/article/10.1086/301102/meta and
 *   https://link.springer.com/article/10.1023%2FA%3A1021149313347 .
 */
template <typename Particles, typename Interaction, typename Regularitor>
class ReguSystem : public particleSystem<Particles, Interaction>
{
public:
    /* Typedef */
    using Base = particleSystem<Particles, Interaction>;
    
    using typename Base::type;
    
    using typename Base::ActiveScalarArray;
    
    using Scalar = typename type::Scalar;
    
    using Vector = typename type::Vector;
    
    using VectorArray = typename type::VectorArray;
    /* Typedef */
    
    using Base::act;
    
    using Base::partc;
    
    /**  @brief Omega interface. Reference to partc.omega*/
    inline Scalar& omega(){ return partc.omega(); }
    
    /**  @brief Bindine energy interface. Reference to partc.bindE*/
    inline Scalar& bindE(){ return partc.bindE(); }

    /**  @brief Advance position one step with current velocity.
     *
     *  Advance position array and physical time one step with current integration step size and velocity.
     *  @param  timeStepSize Integration step size, will be transfered to physical time in the function.
     */
    inline void drift(Scalar stepSize)
    {
        Scalar physicalTime = regular.getPhysicalPosTime(partc, stepSize);
        
        partc.advancePos(partc.vel(), physicalTime);
        partc.advanceTime(physicalTime);
    }
    
    /** @brief Advance velocity one step with current acceleration.
     *
     *  Advance velocity array one step with current integration step size and accelerations.
     *  @param  stepSize Integration step size, will be transfered to physical time in the function.
     */
    inline void kick(Scalar stepSize)
    {
        Scalar physicalTime = regular.getPhysicalVelTime(partc, stepSize);
        
        this->act.zeroTotalAcc();
        
        this->act.calcuVelIndepAcc(partc.mass(), partc.pos(), partc.vel());
        this->act.calcuExtVelIndepAcc(partc.mass(), partc.pos(), partc.vel());
        
        advanceVels<Interaction::isVelDep>(physicalTime);
    }
    
    /** @brief Interface to rescale the time.
     *
     *  Interace used by dynamic system. Transfer integration time to physical time.
     *  @return The phsyical time.
     */
    Scalar timeScale(Scalar scale)
    {
        return regular.getPhysicalPosTime(partc, scale);
    }

private:
    /** @brief Regularization interface.*/
    Regularitor regular;
    
private:
    /** @brief SFINAE version of kick() of velocity independent force */
    template<bool isVelDep>
    inline typename std::enable_if<isVelDep==false>::type
    advanceVels(Scalar stepSize)
    {
        act.calcuTotalAcc();
        partc.advanceOmega(act.velIndepAcc(), partc.vel(), 0.5*stepSize);
        partc.advanceVel(act.totalAcc(), stepSize);
        partc.advanceOmega(act.velIndepAcc(), partc.vel(), 0.5*stepSize);
    }
    
    /** @brief SFINAE version of kick() of velocity dependent force */
    template<bool isVelDep>
    inline typename std::enable_if<isVelDep==true>::type
    advanceVels(Scalar stepSize)
    {
        advanceAuxiVel(stepSize*0.5);
        advanceVel(stepSize);
        partc.advanceOmega(act.velIndepAcc(), partc.auxiVel(), stepSize);
        partc.advanceBindE(act.velDepAcc(),   partc.auxiVel(), stepSize);
        advanceAuxiVel(stepSize*0.5);
    }
    
    /** @brief Advance real velocity one step with current acceleration.
     *
     *  Advance velocity array one step with current integration step size and accelerations.
     *  @param  stepSize Integration step size, will be transfered to physical time in the function.
     */
    inline void advanceVel(Scalar stepSize)
    {
        act.calcuVelDepAcc(partc.mass(), partc.pos(), partc.auxiVel());
        act.calcuExtVelDepAcc(partc.mass(), partc.pos(), partc.auxiVel());
        
        act.calcuTotalAcc();
        partc.advanceVel(act.totalAcc(), stepSize);
    }
    
    /** @brief Advance auxilary velocity one step with current acceleration.
     *
     *  Advance auxi velocity array one step with current integration step size and accelerations.
     *  @param  stepSize Integration step size, will be transfered to physical time in the function.
     */
    inline void advanceAuxiVel(Scalar stepSize)
    {
        act.calcuVelDepAcc(partc.mass(), partc.pos(), partc.vel());
        act.calcuExtVelDepAcc(partc.mass(), partc.pos(), partc.vel());
        
        act.calcuTotalAcc();
        partc.advanceAuxiVel(act.totalAcc(), stepSize);
    }
};

#endif

