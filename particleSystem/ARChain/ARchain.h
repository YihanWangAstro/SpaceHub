
#ifndef ARCHAIN_H
#define ARCHAIN_H
#include "../../particleSystem.h"
#include "../../libs.h"
#include "chain.h"
#include <cstring>
#include <algorithm>

/**  @brief Algorithmatic Regularization chain System.
 *
 *   See details in https://arxiv.org/abs/0709.3367 .
 */
template <typename Particles, typename Interaction, typename Regularitor>
class ARchain : public ReguSystem<Particles, Interaction, Regularitor>
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
    
    using Base::regular;
    /** @brief Advance velocity one step with current acceleration.
     *
     *  Advance velocity array one step with current integration step size and accelerations.
     *  @param  stepSize Integration step size, will be transfered to physical time in the function.
     */
    inline void kick(Scalar stepSize)
    {
        Scalar physicalTime = regular.getPhysicalVelTime(partc, stepSize);
        
        this->act.zeroTotalAcc();
        
        this->act.calcuVelIndepAcc(partc.mass(), partc.pos(), partc.vel(), partc.chainPos(), partc.chainVel());
        this->act.calcuExtVelIndepAcc(partc.mass(), partc.pos(), partc.vel());
        
        this->advanceVels<Interaction::isVelDep>(physicalTime);
    }
    
private:
    /** @brief Advance real velocity one step with current acceleration.
     *
     *  Advance velocity array one step with current integration step size and accelerations.
     *  @param  stepSize Integration step size, will be transfered to physical time in the function.
     */
    inline void advanceRrealVel(Scalar stepSize)
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

/*
template<typename Interaction, typename EvolvedData, typename Regularitor>
std::istream& ARchain<Interaction, EvolvedData, Regularitor>::read(std::istream& input)
{
    input >> this->dynState.time;
    size_t id;

    for(size_t i = 0 ; i < size() ; ++i)
        input >> id >> this->tp[i] >> this->m[i] >> this->rad[i] >> this->dynState.pos[i] >> this->dynState.vel[i];

    memset(&(this->acc[0]), 0, sizeof(Vector)*size());
    this->dynState.moveToCentralMassCoords(this->m);
    this->dynState.initAddiVariable(this->m);
    chain::getChainIndex(this->dynState.pos, chainIndex);
    this->dynState.toChain(chain, chainIndex);
#ifdef KAHAN_SUMMATION
    roundoffErr.setZero();
#endif
    return input;
}*/

/** @brief Load data from a plain array.
 *
 *  Interface usded by integrator and ODE iterator. Load data from a plain array processed by itegrator and
 *  iterator to chain data. Then synchrosize Cartesian data and update the chain. Derived class could overload
 *  this function to additional process.
 *
 *  @param data Plain scalar array.
 */
/*template<typename Interaction, typename EvolvedData, typename Regularitor>
void ARchain<Interaction, EvolvedData, Regularitor>::afterIterProcess()
{
    chain.auxiVel = chain.vel;
    chain.toCartesian(this->dynState, chainIndex);
    this->dynState.moveToCentralMassCoords(this->m);
    updateChain();
#ifdef KAHAN_SUMMATION
    roundoffErr.setZero();
#endif
}
*/


/** @brief Update chain index and chain.
 *
 *  Due to the evolution, the relative position of particles may vary with time. This function update the chain
 *  index and chain if needed.
 */
/*
template<typename Interaction, typename EvolvedData, typename Regularitor>
void ARchain<Interaction, EvolvedData, Regularitor>::updateChain()
{
    IndexArray newIndex;
    chain::getChainIndex(this->dynState.pos, newIndex);

    if(chain::IsDiff(chainIndex, newIndex))
    {
        //Update chain index and update new position chain with old chain
        chain::updateChain(this->chain.pos, chainIndex, newIndex);
        chainIndex = newIndex;
        chain::synChain(this->dynState.vel, chain.vel, chainIndex);
        chain::synChain(this->dynState.auxiVel, chain.auxiVel, chainIndex);
    }
}*/



#endif

