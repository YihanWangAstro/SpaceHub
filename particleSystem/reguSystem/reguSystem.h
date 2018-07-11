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
    
    template<typename T, size_t S>
    using Container   = typename Base::template Container<T, S>;
    
    using Scalar      = typename Base::Scalar;
    
    using Vector      = typename Base::Vector;
    
    using VectorArray = typename Base::VectorArray;
    
    using ScalarArray = typename Base::ScalarArray;
    
    using IntArray    = typename Base::IntArray;
    
    using SizeArray   = typename Base::SizeArray;
    
    using ActiveScalarArray = typename Base::ActiveScalarArray;
    /* Typedef */
    
    /**  @brief Omega interface. Reference to partc.omega*/
    inline Scalar& omega(){ return this->partc.omega(); }
    
    /**  @brief Bindine energy interface. Reference to partc.bindE*/
    inline Scalar& bindE(){ return this->partc.bindE(); }

    /**  @brief Advance position one step with current velocity.
     *
     *  Advance position array and physical time one step with current integration step size and velocity.
     *  @param  timeStepSize Integration step size, will be transfered to physical time in the function.
     */
    inline void drift(Scalar stepSize)
    {
        Scalar physicalTime = regular.getPhysicalPosTime(this->partc, stepSize);
        
        this->partc.advancePos(this->partc.vel(), physicalTime);
        this->partc.advanceTime(physicalTime);
    }
    
    /** @brief Advance velocity one step with current acceleration.
     *
     *  Advance velocity array one step with current integration step size and accelerations.
     *  @param  stepSize Integration step size, will be transfered to physical time in the function.
     */
    inline void kick(Scalar stepSize)
    {
        Scalar physicalTime = regular.getPhysicalVelTime(this->partc, stepSize);
        
        this->act.zeroTotalAcc();
        
        this->act.calcuVelIndepAcc(this->partc.mass(), this->partc.pos(), this->partc.vel());
        this->act.calcuExtVelIndepAcc(this->partc.mass(), this->partc.pos(), this->partc.vel());
        
        advanceVels<Interaction::isVelDep>(physicalTime);
    }
    
    /** @brief Interface to rescale the time.
     *
     *  Interace used by dynamic system. Transfer integration time to physical time.
     *  @return The phsyical time.
     */
    Scalar timeScale(Scalar scale)
    {
        return regular.getPhysicalPosTime(this->partc, scale);
    }

private:
    /** @brief Regularization interface.*/
    Regularitor regular;
    
private:
    template<bool isVelDep>
    inline typename std::enable_if<isVelDep==false>::type
    advanceVels(Scalar stepSize)
    {
        this->act.calcuTotalAcc();
        this->partc.advanceOmega(this->act.velIndepAcc(), this->partc.vel(), 0.5*stepSize);
        this->partc.advanceVel(this->act.totalAcc(), stepSize);
        this->partc.advanceOmega(this->act.velIndepAcc(), this->partc.vel(), 0.5*stepSize);
    }
    
    template<bool isVelDep>
    inline typename std::enable_if<isVelDep==true>::type
    advanceVels(Scalar stepSize)
    {
        advanceAuxiVel(stepSize*0.5);
        advanceVel(stepSize);
        this->partc.advanceOmega(this->act.velIndepAcc(), this->partc.auxiVel(), stepSize);
        this->partc.advanceBindE(this->act.velDepAcc(),   this->partc.auxiVel(), stepSize);
        advanceAuxiVel(stepSize*0.5);
    }
    
    /** @brief Advance real velocity one step with current acceleration.
     *
     *  Advance velocity array one step with current integration step size and accelerations.
     *  @param  stepSize Integration step size, will be transfered to physical time in the function.
     */
    inline void advanceVel(Scalar stepSize)
    {
        this->act.calcuVelDepAcc(this->partc.mass(), this->partc.pos(), this->partc.auxiVel());
        this->act.calcuExtVelDepAcc(this->partc.mass(), this->partc.pos(), this->partc.auxiVel());
        
        this->act.calcuTotalAcc();
        this->partc.advanceVel(this->act.totalAcc(), stepSize);
    }
    
    /** @brief Advance auxilary velocity one step with current acceleration.
     *
     *  Advance auxi velocity array one step with current integration step size and accelerations.
     *  @param  stepSize Integration step size, will be transfered to physical time in the function.
     */
    inline void advanceAuxiVel(Scalar stepSize)
    {
        this->act.calcuVelDepAcc(this->partc.mass(), this->partc.pos(), this->partc.vel());
        this->act.calcuExtVelDepAcc(this->partc.mass(), this->partc.pos(), this->partc.vel());
        
        this->act.calcuTotalAcc();
        this->partc.advanceAuxiVel(this->act.totalAcc(), stepSize);
    }
};

#endif

