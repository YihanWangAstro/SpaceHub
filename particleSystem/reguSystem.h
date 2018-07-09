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
#include "../particleSystem.h"
#include "../interaction/interaction.h"
#include "../libs.h"
#include <cstring>
#include <algorithm>

/**  @brief Regularized particle System.
 *
 *   Regularied particle system. See details in https://link.springer.com/article/10.1023%2FA%3A1008368322547 ,
 *   http://iopscience.iop.org/article/10.1086/301102/meta and
 *   https://link.springer.com/article/10.1023%2FA%3A1021149313347 .
 */
template<typename Interaction, typename Particles, typename Regularitor>
class reguSystem : public particleSystem<Particles>
{
public:
    using Base = particleSystem<Particles>;
    using particleSystem<Particles>::size;
    using particleSystem<Particles>::volume;

    /**  @brief Omega interface. Reference to partc.omega*/
    inline typename Base::Scalar& omega(){ return this->partc.omega(); }
    
    /**  @brief Bindine energy interface. Reference to partc.bindE*/
    inline typename Base::Scalar& bindE(){ return this->partc.bindE(); }
    
    /** @brief Advance position one step with current velocity. */
    void advancePos(typename Base::Scalar timeStepSize);

    /** @brief Advance velocity one step with current acceleration. */
    void advanceVel(typename Base::Scalar timeStepSize);

    /** @brief Interface to rescale the time.
     *  @note  Overload base class timeScale().
     */
    typename Base::Scalar timeScale(typename Base::Scalar scale);
    
    /** @brief After process after iteration*/
    void afterIterProcess();
    
    /** @brief Overload operator << */
    friend std::ostream& operator<<(std::ostream& output, const reguSystem& sys)
    {
        return sys.write(output);
    }
    /** @brief Input from istream */
    friend std::istream& operator>>(std::istream& input, reguSystem& sys)
    {
        return sys.read(input);
    }
    
private:
    /** @brief Velocity independent acceleration array*/
    typename Base::VectorArray velIndepAcc;

    /** @brief Velocity dependent acceleration array*/
    typename Base::VectorArray velDepAcc;

    /** @brief Velocity dependent pair force functor*/
    Interaction velDepForce;

    /** @brief Regularization interface.*/
    Regularitor regular;

private:
    /** @brief Advance velocity with current acceleration.*/
    void kickVel(typename Base::Scalar stepSize);

    /** @brief Advance auxiliar velocity with current acceleration.*/
    void kickAuxiVel(typename Base::Scalar stepSize);

    /** @brief Update velocity dependent acceleration with given velocity. Then update the total acceleration.*/
    void updateAccWith(typename Base::VectorArray& vel);

    /** @brief Update velocity independent acceleration.*/
    void updateVelIndepAcc();
    
};

/**  @brief Advance position one step with current velocity.
  *
  *  Advance position array and physical time one step with current integration step size and velocity.
  *  @param  timeStepSize Integration step size, will be transfered to physical time in the function.
 */
template<typename Interaction, typename Particles, typename Regularitor>
void reguSystem<Interaction, Particles, Regularitor>::advancePos(typename Base::Scalar timeStepSize)
{
    Scalar physicalTime = regular.getPhysicalPosTime(this->partc, timeStepSize);
    this->partc.advancePos(physicalTime);
    this->partc.advanceTime(physicalTime);
}

/** @brief Advance velocity one step with current acceleration.
 *
 *  Advance velocity and auxiliary velocity array one step with current integration step size and accelerations.
 *  @param  timeStepSize Integration step size, will be transfered to physical time in the function.
 */
template<typename Interaction, typename Particles, typename Regularitor>
void reguSystem<Interaction, Particles, Regularitor>::advanceVel(typename Base::Scalar timeStepSize)
{
    Scalar physicalTime = regular.getPhysicalVelTime(this->partc, timeStepSize);
    Scalar halfTime = 0.5 * physicalTime;
    updateVelIndepAcc();
    
    updateAccWith(this->dynState.vel);
    this->partc.advanceAuxiVel(this->acc, halfTime);
    updateAccWith(this->dynState.auxiVel);
    this->partc.advanceBindE(,physicalTime);
    this->partc.advanceOmega(,physicalTime);
    this->partc.advanceAuxiVel(this->acc, physicalTime);
    updateAccWith(this->dynState.vel);
    this->partc.advanceAuxiVel(this->acc, halfTime);

}

/** @brief Process data after one step iteration.
 *
 *  Interface usded by ODE iterator. Synchronize auxiliary velocity and set round off error to zero
 *  if Kahan summation is used.
 */
template<typename Interaction, typename Particles, typename Regularitor>
void reguSystem<Interaction, Particles, Regularitor>::afterIterProcess()
{
#ifdef KAHAN_SUMMATION
    this->partc.resetErr();
#endif
}

/** @brief Interface to rescale the time.
 *
 *  Interace used by dynamic system. Transfer integration time to physical time.
 *  @return The phsyical time.
 */
template<typename Interaction, typename Particles, typename Regularitor>
typename Base::Scalar reguSystem<Interaction, Particles, Regularitor>::timeScale(typename Base::Scalar scale)
{
    return regular.getPhysicalPosTime(this->partc, scale);
}


/** @brief Advance regularization variable omega.
 *
 *  Advance omega with velocity independent acceleration and auxiliar velocity with physical time step size.
 *  @param stepSize Physical time step.
 */
template<typename Interaction, typename Particles, typename Regularitor>
void reguSystem<Interaction, Particles, Regularitor>::advanceOmega(typename Base::Scalar stepSize)
{
    Scalar dOmega = 0;

    for(size_t i = 0 ; i < size() ; ++i)
        dOmega += (this->velIndepAcc[i] * this->dynState.auxiVel[i]) * (this->m[i]);

    this->partc.advanceOmega( dOmega * stepSize);
}

/** @brief Advance regularization variable bindE.
 *
 *  Advance bindE with velocity dependent acceleration and auxiliar velocity with physical time step size.
 *  @param stepSize Physical time step.
 */
template<typename Interaction, typename Particles, typename Regularitor>
void reguSystem<Interaction, Particles, Regularitor>::advanceB(typename Base::Scalar stepSize)
{
    Scalar dE = 0;

    for(size_t i = 0 ; i < size() ; ++i)
        dE -= (this->velDepAcc[i] * this->dynState.auxiVel[i]) * (this->m[i]);

    this->partc.advanceBindE(dE * stepSize);
}

/** @brief Update velocity independent acceleration.
 *
 *  Update velocity independent acceleration 'velIndepAcc' with Newtonian interaction.
 */
template<typename Interaction, typename Particles, typename Regularitor>
void reguSystem<Interaction, Particles, Regularitor>::updateVelIndepAcc()
{
    Vector dr(0.0, 0.0, 0.0);
    Scalar inv_r  = 1;
    Scalar inv_r3 = 1;
    memset(&(velIndepAcc[0]), 0, sizeof(Vector)*size());

    for(size_t i = 0 ; i < size() ; ++i)
    {
        for(size_t j = i + 1 ; j < size() ; ++j)
        {
            dr     = this->dynState.pos[j] - this->dynState.pos[i];
            inv_r  = dr.reNorm();
            inv_r3 = inv_r * inv_r * inv_r;
            velIndepAcc[i] += dr * (inv_r3 * this->m[j]);
            velIndepAcc[j] -= dr * (inv_r3 * this->m[i]);
        }
    }
}

/** @brief Update velocity dependent acceleration with given velocity. Then update the total acceleration.
 *
 *  Update the velocity dependent accelaration variable 'velDepAcc' with given velocity and velocity dependent pair
 *  force 'velDepForce'. Then update the total acceleration 'acc' by adding 'velIndepAcc' and 'velDepAcc'.
 */
template<typename Interaction, typename Particles, typename Regularitor>
void reguSystem<Interaction, Particles, Regularitor>::updateAccWith(typename Base::VectorArray& velocity)
{
    memset(&(velDepAcc[0]), 0, sizeof(Vector)*size());
    Vector dr;
    Vector dv;

    for(size_t i = 0 ; i < size() ; ++i)
    {
        for(size_t j = i + 1 ; j < size() ; ++j)
        {
            dr = this->dynState.pos[i] - this->dynState.pos[j];
            dv = velocity[i] - velocity[j];
            velDepForce(this->m[i], this->m[j], dr, dv,
                        velocity[i], velocity[j],
                        velDepAcc[i], velDepAcc[j]);
        }
    }

    for(size_t i = 0 ; i < size() ; ++i)
        this->acc[i] = velIndepAcc[i] + velDepAcc[i];
}

/*==============================Specialization================================*/
/**  @brief Partial specialization of reguSystem for velocity independent system. */
template<typename Particles, typename Regularitor>
class reguSystem<interact::Newtonian<typename Particles::Scalar>, Particles, Regularitor> : public particleSystem<Particles>
{
public:
    
    using particleSystem<Particles>::size;
    using particleSystem<Particles>::volume;
    
    /**  @brief Omega interface. Reference to dynState.omega*/
    inline Scalar& omega(){ return this->dynState.omega; }
    
    /**  @brief Bindine energy interface. Reference to dynState.bindE*/
    inline Scalar& bindE(){ return this->dynState.bindE; }
    
    /** @brief Advance position one step with current velocity. */
    void advancePos(Scalar timeStepSize);
    
     /** @brief Advance velocity one step with current acceleration. */
    void advanceVel(Scalar timeStepSize);
    
    /** @brief Interface to rescale the time.
     *  @note  Overload base class timeScale().
     */
    Scalar timeScale(Scalar scale);
    
    /** @brief After process after iteration*/
    void afterIterProcess();
    
    /** @brief Overload operator << */
    friend std::ostream& operator<<(std::ostream& output, const reguSystem& sys)
    {
        return sys.write(output);
    }
    /** @brief Input from istream */
    friend std::istream& operator>>(std::istream& input, reguSystem& sys)
    {
        return sys.read(input);
    }
    
private:
#ifdef KAHAN_SUMMATION
    /**  @brief Co-evolved Data used by Kahan summation.
     *
     *   See details in https://en.wikipedia.org/wiki/Kahan_summation_algorithm .
     */
    Particles  roundoffErr;
#endif
    
    /** @brief Regularization interface.*/
    Regularitor regular;
    
private:
    
    /** @brief Advance regularization variable omega.*/
    void advanceOmega(Scalar stepSize);
    
    /** @brief Update velocity independent acceleration.*/
    void updateVelIndepAcc();
    
};

/**  @brief Advance position one step with current velocity.
 *
 *  Advance position array and physical time one step with current integration step size and velocity.
 *  @param  timeStepSize Integration step size, will be transfered to physical time in the function.
 */
template<typename Particles, typename Regularitor>
void reguSystem<interact::Newtonian<typename Particles::Scalar>, Particles, Regularitor>::advancePos(Scalar timeStepSize)
{
    Scalar physicalTime = regular.getPhysicalPosTime(this->m, this->dynState, timeStepSize);
#ifdef KAHAN_SUMMATION
    KahanAdvance(this->dynState.pos, this->dynState.vel, roundoffErr.pos, physicalTime);
#else
    advanceVariable(this->dynState.pos, this->dynState.vel, physicalTime);
#endif
    this->dynState.time += physicalTime;
}

/** @brief Advance velocity one step with current acceleration.
 *
 *  Advance velocity array one step with current integration step size and accelerations.
 *  @param  timeStepSize Integration step size, will be transfered to physical time in the function.
 */
template<typename Particles, typename Regularitor>
void reguSystem<interact::Newtonian<typename Particles::Scalar>, Particles, Regularitor>::advanceVel(Scalar timeStepSize)
{
    Scalar physicalTime = regular.getPhysicalVelTime(this->m, this->dynState, timeStepSize);
    updateVelIndepAcc();
    advanceOmega(0.5 * physicalTime);
#ifdef KAHAN_SUMMATION
    KahanAdvance(this->dynState.vel, this->acc, roundoffErr.vel, physicalTime);
#else
    advanceVariable(this->dynState.vel, this->acc, physicalTime);
#endif
    advanceOmega(0.5 * physicalTime);
}

/** @brief Process data after one step iteration.
 *
 *  Interface usded by ODE iterator. Synchronize auxiliary velocity and set round off error to zero
 *  if Kahan summation is used.
 */
template<typename Particles, typename Regularitor>
void reguSystem<interact::Newtonian<typename Particles::Scalar>, Particles, Regularitor>::afterIterProcess()
{
#ifdef KAHAN_SUMMATION
    roundoffErr.setZero();
#endif
}

/** @brief Interface to rescale the time.
 *
 *  Interace used by dynamic system. Transfer integration time to physical time.
 *  @return The phsyical time.
 */
template<typename Particles, typename Regularitor>
typename Particles::Scalar reguSystem<interact::Newtonian<typename Particles::Scalar>, Particles, Regularitor>::timeScale(
    Scalar scale)
{
    return regular.getPhysicalPosTime(this->m, this->dynState, scale);
}

/** @brief Advance regularization variable omega.
 *
 *  Advance omega with velocity independent acceleration and auxiliar velocity with physical time step size.
 *  @param stepSize Physical time step.
 */
template<typename Particles, typename Regularitor>
void reguSystem<interact::Newtonian<typename Particles::Scalar>, Particles, Regularitor>::advanceOmega(Scalar stepSize)
{
    Scalar dOmega = 0;

    for(size_t i = 0 ; i < size() ; ++i)
        dOmega += (this->acc[i] * this->dynState.vel[i]) * (this->m[i]);

#ifdef KAHAN_SUMMATION
    KahanAdvance(this->dynState.omega, dOmega * stepSize, roundoffErr.omega);
#else
    this->dynState.omega += dOmega * stepSize;
#endif
}

/** @brief Update velocity independent acceleration.
 *
 *  Update velocity independent acceleration 'acc' with Newtonian interaction.
 */
template<typename Particles, typename Regularitor>
void reguSystem<interact::Newtonian<typename Particles::Scalar>, Particles, Regularitor>::updateVelIndepAcc()
{
    Vector dr(0.0, 0.0, 0.0);
    Scalar inv_r  = 1;
    Scalar inv_r3 = 1;
    std::for_each(this->acc.begin(), this->acc.end(), [](vec3<Scalar>& v)
    {
        v.setZero();
    });

    for(size_t i = 0 ; i < size() ; ++i)
    {
        for(size_t j = i + 1 ; j < size() ; ++j)
        {
            dr     = this->dynState.pos[j] - this->dynState.pos[i];
            inv_r  = dr.reNorm();
            inv_r3 = inv_r * inv_r * inv_r;
            this->acc[i] += dr * (inv_r3 * this->m[j]);
            this->acc[j] -= dr * (inv_r3 * this->m[i]);
        }
    }
}
#endif

