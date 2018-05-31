
#ifndef ARCHAIN_H
#define ARCHAIN_H
#include "../particleSystem.h"
#include "../libs.h"
#include "chain.h"
#include <cstring>
#include <algorithm>

/**  @brief Algorithmatic Regularization chain System.
 *
 *   See details in https://arxiv.org/abs/0709.3367 .
 */
template<typename Interaction, typename EvolvedData, typename Regularitor>
class ARchain : public particleSystem<EvolvedData>
{
public:
    
    typedef typename EvolvedData::Scalar      Scalar;
    typedef typename EvolvedData::Vector      Vector;
    typedef typename EvolvedData::VectorArray VectorArray;
    typedef typename EvolvedData::ScalarArray ScalarArray;
    typedef typename EvolvedData::PlainArray  PlainArray;
    typedef typename EvolvedData::IndexArray  IndexArray;
    using particleSystem<EvolvedData>::size;
    using particleSystem<EvolvedData>::volume;
    
    /**  @brief Omega interface. Reference to dynState.omega*/
    inline Scalar& omega(){ return this->dynState.omega; }
    
    /**  @brief Bindine energy interface. Reference to dynState.bindE*/
    inline Scalar& bindE(){ return this->dynState.bindE; }
    
    /** @brief Advance position one step with current velocity. */
    void advancePos(Scalar timeStepSize);
    
    /** @brief Advance velocity one step with current acceleration. */
    void advanceVel(Scalar timeStepSize);
    
    /** @brief Input data from standard c++ istream.
     *  @note  Overload base class read().
     */
    std::istream& read(std::istream&);
    
    /** @brief Interface to rescale the time.
     *  @note  Overload base class timeScale().
     */
    Scalar timeScale(Scalar scale);
    
    /** @brief Transfer evolved data to a plain array.
     *  @param arr The destination plain array.
     */
    void flatten(PlainArray& arr)
    {
        chain.flatten(arr);
    }
    
    /** @brief Load data from a plain array.
     *  @param arr The plain array data.
     */
    void loadFlatten(PlainArray& arr)
    {
        chain.loadFlatten(arr);
    }
    
    /** @brief After process after iteration*/
    void afterIterProcess();
    
    /** @brief Overload operator << */
    friend std::ostream& operator<<(std::ostream& output, const ARchain& sys)
    {
        return sys.write(output);
    }
    /** @brief Input from istream */
    friend std::istream& operator>>(std::istream& input, ARchain& sys)
    {
        return sys.read(input);
    }
private:
#ifdef KAHAN_SUMMATION
    /**  @brief Co-evolved Data used by Kahan summation.
     *
     *   See details in https://en.wikipedia.org/wiki/Kahan_summation_algorithm .
     */
    EvolvedData roundoffErr;
#endif
    /**  @brief Evolved variables in chain coordinates*/
    EvolvedData chain;
    
    /**  @brief Mapping index from Cartesian coordinates to chain coordinates*/
    IndexArray  chainIndex;
    
    /** @brief Velocity independent acceleration array in Cartesian coordiantes*/
    VectorArray velIndepAcc;
    
    /** @brief Velocity independent acceleration array in chain coordiantes*/
    VectorArray chainVelIndepAcc;
    
    /** @brief Velocity dependent acceleration array in Cartesian coordiantes*/
    VectorArray velDepAcc;
    
    /** @brief Velocity dependent acceleration array in chain coordiantes*/
    VectorArray chainVelDepAcc;
    
    /** @brief Velocity dependent pair force functor*/
    Interaction velDepForce;
    
    /** @brief Regularization interface.*/
    Regularitor regular;
    
    /** @brief Advance regularization variable omega.*/
    void advanceOmega(Scalar stepSize);
    
    /** @brief Advance regularization variable bindE.*/
    void advanceB(Scalar stepSize);
    
    /** @brief Advance velocity with current acceleration.*/
    void kickVel(Scalar stepSize);
    
    /** @brief Advance auxiliary velocity with current acceleration.*/
    void kickAuxiVel(Scalar stepSize);
    
    /** @brief Update velocity dependent acceleration with given velocity. Then update the total acceleration.*/
    void updateAccWith(VectorArray& vel, VectorArray& chainVel);
    
    /** @brief Update velocity independent acceleration.*/
    void updateVelIndepAcc();
    
    /** @brief Update the chain data based on current Cartesian coordinates.*/
    void updateChain();
    
};

/** @brief Advance chain position one step with current velocity.
 *
 *  Advance chain position array and physical time one step with current integration step size and velocity.
 *  @param  timeStepSize Integration step size, will be transfered to physical time in the function.
 */
template<typename Interaction, typename EvolvedData, typename Regularitor>
void ARchain<Interaction, EvolvedData, Regularitor>::advancePos(Scalar timeStepSize)
{
    Scalar physicalTime = regular.getPhysicalPosTime(this->m, this->dynState, timeStepSize);
#ifdef KAHAN_SUMMATION
    KahanAdvance(chain.pos, chain.vel, roundoffErr.pos, physicalTime);
#else
    advanceVariable(chain.pos, chain.vel, physicalTime);
#endif
    chain::synCartesian(chain.pos, this->dynState.pos, chainIndex);
    MoveToCentralMassCoordinate(this->m, this->dynState.pos);
    this->chain.time += physicalTime;
    this->dynState.time = this->chain.time;
}

/** @brief Advance chain velocity one step with current acceleration.
 *
 *  Advance chain velocity and chian auxiliary velocity array one step with current integration step size
 *  and accelerations.
 *  @param  timeStepSize Integration step size, will be transfered to physical time in the function.
 */
template<typename Interaction, typename EvolvedData, typename Regularitor>
void ARchain<Interaction, EvolvedData, Regularitor>::advanceVel(Scalar timeStepSize)
{
    Scalar physicalTime = regular.getPhysicalVelTime(this->m, this->dynState, timeStepSize);
    Scalar halfTime     = 0.5 * physicalTime;
    updateVelIndepAcc();
    updateAccWith(this->dynState.vel, chain.vel);
    kickAuxiVel(halfTime);
    updateAccWith(this->dynState.auxiVel, chain.auxiVel);
    advanceB(physicalTime);
    advanceOmega(physicalTime);
    kickVel(physicalTime);
    updateAccWith(this->dynState.vel, chain.vel);
    kickAuxiVel(halfTime);
}

/** @brief Input data from standard c++ istream.
 *
 *  Implement of CRTP '>>' method.
 */
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
}

/** @brief Load data from a plain array.
 *
 *  Interface usded by integrator and ODE iterator. Load data from a plain array processed by itegrator and
 *  iterator to chain data. Then synchrosize Cartesian data and update the chain. Derived class could overload
 *  this function to additional process.
 *
 *  @param data Plain scalar array.
 */
template<typename Interaction, typename EvolvedData, typename Regularitor>
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

/** @brief Interface to rescale the time.
 *
 *  Interace used by dynamic system. Transfer integration time to physical time.
 *  @return The phsyical time.
 */
template<typename Interaction, typename EvolvedData, typename Regularitor>
typename EvolvedData::Scalar ARchain<Interaction, EvolvedData, Regularitor>::timeScale(Scalar scale)
{
    return regular.getPhysicalPosTime(this->m, this->dynState, scale);
}

/** @brief Advance chain velocity with current acceleration.
 *
 *  Advance the chain velocity with current total acceleration variable 'acc' and synchronize the Cartesian data.
 *  @param stepSize Integration step size.
 */
template<typename Interaction, typename EvolvedData, typename Regularitor>
void ARchain<Interaction, EvolvedData, Regularitor>::kickVel(Scalar stepSize)
{
#ifdef KAHAN_SUMMATION
    KahanAdvance(chain.vel, this->acc, roundoffErr.vel, stepSize);
#else
    advanceVariable(chain.vel, this->acc, stepSize);
#endif
    chain::synCartesian(chain.vel, this->dynState.vel, chainIndex);
    MoveToCentralMassCoordinate(this->m, this->dynState.vel);
}

/** @brief Advance chain auxiliary velocity with current acceleration.
 *
 *  Advance the chain auxiliary velocity with current total acceleration variable 'acc' and synchronize the
 *  Cartesian data.
 *  @param stepSize Integration step size.
 */
template<typename Interaction, typename EvolvedData, typename Regularitor>
void ARchain<Interaction, EvolvedData, Regularitor>::kickAuxiVel(Scalar stepSize)
{
#ifdef KAHAN_SUMMATION
    KahanAdvance(chain.auxiVel, this->acc, roundoffErr.auxiVel, stepSize);
#else
    advanceVariable(chain.auxiVel, this->acc, stepSize);
#endif
    chain::synCartesian(chain.auxiVel, this->dynState.auxiVel, chainIndex);
    MoveToCentralMassCoordinate(this->m, this->dynState.auxiVel);
}

/** @brief Advance regularization variable omega.
 *
 *  Advance omega with velocity independent acceleration and auxiliar velocity with physical time step size.
 *  @param stepSize Physical time step.
 *  @note  Update the omega in chain, which will be processed by the integrator.
 */
template<typename Interaction, typename EvolvedData, typename Regularitor>
void ARchain<Interaction, EvolvedData, Regularitor>::advanceOmega(Scalar stepSize)
{
    Scalar dOmega = 0;

    for(size_t i = 0 ; i < size() ; ++i)
        dOmega += (this->velIndepAcc[i] * this->dynState.auxiVel[i]) * (this->m[i]);

#ifdef KAHAN_SUMMATION
    KahanAdvance(this->chain.omega, dOmega * stepSize, roundoffErr.omega);
#else
    this->chain.omega += dOmega * stepSize;
#endif
    this->dynState.omega = this->chain.omega;
}

/** @brief Advance regularization variable bindE.
 *
 *  Advance bindE with velocity dependent acceleration and auxiliar velocity with physical time step size.
 *  @param stepSize Physical time step.
 *  @note  Update the bindE in chain, which will be processed by the integrator.
 */
template<typename Interaction, typename EvolvedData, typename Regularitor>
void ARchain<Interaction, EvolvedData, Regularitor>::advanceB(Scalar stepSize)
{
    Scalar dE = 0;

    for(size_t i = 0 ; i < size() ; ++i)
        dE -= (this->velDepAcc[i] * this->dynState.auxiVel[i]) * (this->m[i]);

#ifdef KAHAN_SUMMATION
    KahanAdvance(this->chain.bindE, dE * stepSize, roundoffErr.bindE);
#else
    this->chain.bindE += dE * stepSize;
#endif
    this->dynState.bindE = this->chain.bindE;
}

/** @brief Update chain index and chain.
 *
 *  Due to the evolution, the relative position of particles may vary with time. This function update the chain
 *  index and chain if needed.
 */
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
}

/** @brief Update velocity independent acceleration.
 *
 *  Update velocity independent acceleration 'velIndepAcc' based on Newtonian interaction with chain data.
 */
template<typename Interaction, typename EvolvedData, typename Regularitor>
void ARchain<Interaction, EvolvedData, Regularitor>::updateVelIndepAcc()
{
    Vector dr(0.0, 0.0, 0.0);
    Scalar inv_r  = 1;
    Scalar inv_r3 = 1;
    memset(&(velIndepAcc[0]), 0, sizeof(Vector)*size());

    for(size_t i = 0 ; i < size() - 1; ++i)
    {
        dr = chain.pos[i];
        inv_r  = dr.reNorm();
        inv_r3 = (inv_r * inv_r * inv_r);
        velIndepAcc[chainIndex[i]]     += dr * (inv_r3 * this->m[chainIndex[i + 1]]);
        velIndepAcc[chainIndex[i + 1]] -= dr * (inv_r3 * this->m[chainIndex[i]]);
    }

    for(size_t i = 0 ; i < size() - 2; ++i)
    {
        dr = chain.pos[i] + chain.pos[i + 1];
        inv_r = dr.reNorm();
        inv_r3 = (inv_r * inv_r * inv_r);
        velIndepAcc[chainIndex[i]]     += dr * (inv_r3 * this->m[chainIndex[i + 2]]);
        velIndepAcc[chainIndex[i + 2]] -= dr * (inv_r3 * this->m[chainIndex[i]]);
    }

    for(size_t i = 0 ; i < size() ; ++i)
    {
        for(size_t j = i + 3 ; j < size(); ++j)
        {
            dr = this->dynState.pos[chainIndex[j]] - this->dynState.pos[chainIndex[i]];
            inv_r = dr.reNorm();
            inv_r3 = (inv_r * inv_r * inv_r);
            velIndepAcc[chainIndex[i]] += dr * (inv_r3 * this->m[chainIndex[j]]);
            velIndepAcc[chainIndex[j]] -= dr * (inv_r3 * this->m[chainIndex[i]]);
        }
    }

    for(size_t i = 0 ; i < size() - 1; ++i)
        chainVelIndepAcc[i] = velIndepAcc[chainIndex[i + 1]] - velIndepAcc[chainIndex[i]];

}

/** @brief Update velocity dependent acceleration with given velocity. Then update the total acceleration.
 *
 *  Update the velocity dependent accelaration variable 'velDepAcc' with given velocity and velocity dependent pair
 *  force 'velDepForce'. Then update the total acceleration 'acc' by adding 'velIndepAcc' and 'velDepAcc'.
 */
template<typename Interaction, typename EvolvedData, typename Regularitor>
void ARchain<Interaction, EvolvedData, Regularitor>::updateAccWith(VectorArray& velocity, VectorArray& chainVel)
{
    memset(&(velDepAcc[0]), 0, sizeof(Vector)*size());
    Vector dr;
    Vector dv;

    for(int i = 1 ; i < size(); ++i)//update chain data(near pair)
    {
        dr = chain.pos[i - 1];
        dv = chainVel[i - 1];
        velDepForce(this->m[chainIndex[i]], this->m[chainIndex[i - 1]], dr, dv,
                    velocity[chainIndex[i]], velocity[chainIndex[i - 1]],
                    velDepAcc[chainIndex[i]], velDepAcc[chainIndex[i - 1]]);
    }

    for(int i = 2 ; i < size(); ++i)//update chain data(second near pair)
    {
        dr = chain.pos[i - 2] + chain.pos[i - 1];
        dv = chainVel[i - 2] + chainVel[i - 1];
        velDepForce(this->m[chainIndex[i]], this->m[chainIndex[i - 2]], dr, dv,
                    velocity[chainIndex[i]], velocity[chainIndex[i - 2]],
                    velDepAcc[chainIndex[i]], velDepAcc[chainIndex[i - 2]]);
    }

    for(int i = 0 ; i < size(); ++i)//update non-chain data
    {
        for(int j = 0 ; j < i - 2; ++j)
        {
            dr = this->dynState.pos[chainIndex[i]] - this->dynState.pos[chainIndex[j]];
            dv = velocity[chainIndex[i]] - velocity[chainIndex[j]];
            velDepForce(this->m[chainIndex[i]], this->m[chainIndex[j]], dr, dv,
                        velocity[chainIndex[i]], velocity[chainIndex[j]],
                        velDepAcc[chainIndex[i]], velDepAcc[chainIndex[j]]);
        }
    }

    for(size_t i = 0 ; i < size() - 1; ++i)
        chainVelDepAcc[i] = velDepAcc[chainIndex[i + 1]] - velDepAcc[chainIndex[i]];

    for(size_t i = 0 ; i < size() ; ++i)
        this->acc[i] = chainVelIndepAcc[i] + chainVelDepAcc[i];
}

/*==============================Specialization================================*/

/** @brief Partial specilization of template ARchain*/
template<typename EvolvedData, typename Regularitor>
class ARchain<interact::Newtonian<typename EvolvedData::Scalar>, EvolvedData, Regularitor> :
    public particleSystem<EvolvedData>
{
public:
    
    typedef typename EvolvedData::Scalar      Scalar;
    typedef typename EvolvedData::Vector      Vector;
    typedef typename EvolvedData::VectorArray VectorArray;
    typedef typename EvolvedData::ScalarArray ScalarArray;
    typedef typename EvolvedData::PlainArray  PlainArray;
    typedef typename EvolvedData::IndexArray  IndexArray;
    using particleSystem<EvolvedData>::size;
    using particleSystem<EvolvedData>::volume;
    
    /**  @brief Omega interface. Reference to dynState.omega*/
    inline Scalar& omega(){ return this->dynState.omega; }
    
    /**  @brief Bindine energy interface. Reference to dynState.bindE*/
    inline Scalar& bindE(){ return this->dynState.bindE; }
    
    /** @brief Advance position one step with current velocity. */
    void advancePos(Scalar timeStepSize);
    
    /** @brief Advance velocity one step with current acceleration. */
    void advanceVel(Scalar timeStepSize);
    
    /** @brief Input data from standard c++ istream.
     *  @note  Overload base class read().
     */
    std::istream& read(std::istream&);
    
    /** @brief Interface to rescale the time.
     *  @note  Overload base class timeScale().
     */
    Scalar timeScale(Scalar scale);
    
    /** @brief Get the number of the particles.
     *  @return The particle number.
     *  @note  Overload base class size().
     */
    constexpr static size_t size()
    {
        return EvolvedData::size();
    }
    
    /** @brief Transfer evolved data to a plain array.
     *  @param arr The destination plain array.
     */
    void flatten(PlainArray& arr)
    {
        chain.flatten(arr);
    }
    
    /** @brief Load data from a plain array.
     *  @param arr The plain array data.
     */
    void loadFlatten(PlainArray& arr)
    {
        chain.loadFlatten(arr);
    }
    
    /** @brief After process after iteration*/
    void afterIterProcess();
 
    /** @brief Overload operator << */
    friend std::ostream& operator<<(std::ostream& output, const ARchain& sys)
    {
        return sys.write(output);
    }
    /** @brief Input from istream */
    friend std::istream& operator>>(std::istream& input, ARchain& sys)
    {
        return sys.read(input);
    }
private:
#ifdef KAHAN_SUMMATION
    /**  @brief Co-evolved Data used by Kahan summation.
     *
     *   See details in https://en.wikipedia.org/wiki/Kahan_summation_algorithm .
     */
    EvolvedData roundoffErr;
#endif
    
    /**  @brief Evolved variables in chain coordinates*/
    EvolvedData chain;
    
    /**  @brief Mapping index from Cartesian coordinates to chain coordinates*/
    IndexArray  chainIndex;
    
    /** @brief Velocity independent acceleration array in Cartesian coordiantes*/
    VectorArray velIndepAcc;
    
    /** @brief Regularization interface.*/
    Regularitor regular;

    /** @brief Advance regularization variable omega.*/
    void advanceOmega(Scalar stepSize);
    
    /** @brief Update velocity independent acceleration.*/
    void updateVelIndepAcc();
    
    /** @brief Update the chain data based on current Cartesian coordinates.*/
    void updateChain();
    
};

/** @brief Advance chain position one step with current velocity.
 *
 *  Advance chain position array and physical time one step with current integration step size and velocity.
 *  @param  timeStepSize Integration step size, will be transfered to physical time in the function.
 */
template<typename EvolvedData, typename Regularitor>
void ARchain<interact::Newtonian<typename EvolvedData::Scalar>, EvolvedData, Regularitor>::advancePos(Scalar timeStepSize)
{
    Scalar physicalTime = regular.getPhysicalPosTime(this->m, this->dynState, timeStepSize);
#ifdef KAHAN_SUMMATION
    KahanAdvance(chain.pos, chain.vel, roundoffErr.pos, physicalTime);
#else
    advanceVariable(chain.pos, chain.vel, physicalTime);
#endif
    chain::synCartesian(chain.pos, this->dynState.pos, chainIndex);
    MoveToCentralMassCoordinate(this->m, this->dynState.pos);
    this->chain.time += physicalTime;
    this->dynState.time = this->chain.time;
}

/** @brief Advance chain velocity one step with current acceleration.
 *
 *  Advance chain velocity one step with current integration step size and accelerations.
 *  @param  timeStepSize Integration step size, will be transfered to physical time in the function.
 */
template<typename EvolvedData, typename Regularitor>
void ARchain<interact::Newtonian<typename EvolvedData::Scalar>, EvolvedData, Regularitor>::advanceVel(Scalar timeStepSize)
{
    Scalar physicalTime = regular.getPhysicalVelTime(this->m, this->dynState, timeStepSize);
    Scalar halfTime     = 0.5 * physicalTime;
    updateVelIndepAcc();
    advanceOmega(halfTime);
#ifdef KAHAN_SUMMATION
    KahanAdvance(chain.vel, this->acc, roundoffErr.vel, physicalTime);
#else
    advanceVariable(chain.vel, this->acc, physicalTime);
#endif
    chain::synCartesian(chain.vel, this->dynState.vel, chainIndex);
    MoveToCentralMassCoordinate(this->m, this->dynState.vel);
    advanceOmega(halfTime);
}

/** @brief Input data from standard c++ istream.
 *
 *  Implement of CRTP '>>' method.
 */
template<typename EvolvedData, typename Regularitor>
std::istream& ARchain<interact::Newtonian<typename EvolvedData::Scalar>, EvolvedData, Regularitor>::read(std::istream& input)
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
}

/** @brief Load data from a plain array.
 *
 *  Interface usded by integrator and ODE iterator. Load data from a plain array processed by itegrator and
 *  iterator to chain data. Then synchrosize Cartesian data and update the chain. Derived class could overload
 *  this function to additional process.
 *
 *  @param data Plain scalar array.
 */
template<typename EvolvedData, typename Regularitor>
void ARchain<interact::Newtonian<typename EvolvedData::Scalar>, EvolvedData, Regularitor>::afterIterProcess()
{
    chain.toCartesian(this->dynState, chainIndex);
    this->dynState.moveToCentralMassCoords(this->m);
    updateChain();
#ifdef KAHAN_SUMMATION
    roundoffErr.setZero();
#endif
}

/** @brief Interface to rescale the time.
 *
 *  Interace used by dynamic system. Transfer integration time to physical time.
 *  @return The phsyical time.
 */
template<typename EvolvedData, typename Regularitor>
typename EvolvedData::Scalar ARchain<interact::Newtonian<typename EvolvedData::Scalar>, EvolvedData, Regularitor>::timeScale(
    Scalar scale)
{
    return regular.getPhysicalPosTime(this->m, this->dynState, scale);
}

/** @brief Advance regularization variable omega.
 *
 *  Advance omega with velocity independent acceleration and velocity with physical time step size.
 *  @param stepSize Physical time step.
 *  @note  Update the omega in chain, which will be processed by the integrator.
 */
template<typename EvolvedData, typename Regularitor>
void ARchain<interact::Newtonian<typename EvolvedData::Scalar>, EvolvedData, Regularitor>::advanceOmega(Scalar stepSize)
{
    Scalar dOmega = 0;

    for(size_t i = 0 ; i < size() ; ++i)
        dOmega += (this->velIndepAcc[i] * this->dynState.vel[i]) * (this->m[i]);

#ifdef KAHAN_SUMMATION
    KahanAdvance(this->chain.omega, dOmega * stepSize, roundoffErr.omega);
#else
    this->chain.omega += dOmega * stepSize;
#endif
    this->dynState.omega = this->chain.omega;
}

/** @brief Update chain index and chain.
 *
 *  Due to the evolution, the relative position of particles may vary with time. This function update the chain
 *  index and chain if needed.
 */
template<typename EvolvedData, typename Regularitor>
void ARchain<interact::Newtonian<typename EvolvedData::Scalar>, EvolvedData, Regularitor>::updateChain()
{
    IndexArray newIndex;
    chain::getChainIndex(this->dynState.pos, newIndex);

    if(chain::IsDiff(chainIndex, newIndex))
    {
        chain::updateChain(this->chain.pos, chainIndex,
                           newIndex);//Update chain index and update new position chain with old chain
        chainIndex = newIndex;
        chain::synChain(this->dynState.vel, chain.vel, chainIndex);
    }
}

/** @brief Update velocity independent acceleration.
 *
 *  Update velocity independent acceleration 'velIndepAcc' based on Newtonian interaction with chain data.
 *  Then calculate the corresponding chain acceleration.
 */
template<typename EvolvedData, typename Regularitor>
void ARchain<interact::Newtonian<typename EvolvedData::Scalar>, EvolvedData, Regularitor>::updateVelIndepAcc()
{
    Vector dr(0.0, 0.0, 0.0);
    Scalar inv_r  = 1;
    Scalar inv_r3 = 1;
    memset(&(velIndepAcc[0]), 0, sizeof(Vector)*size());

    for(size_t i = 0 ; i < size() - 1; ++i)
    {
        dr = chain.pos[i];
        inv_r  = dr.reNorm();
        inv_r3 = (inv_r * inv_r * inv_r);
        velIndepAcc[chainIndex[i]]     += dr * (inv_r3 * this->m[chainIndex[i + 1]]);
        velIndepAcc[chainIndex[i + 1]] -= dr * (inv_r3 * this->m[chainIndex[i]]);
    }

    for(size_t i = 0 ; i < size() - 2; ++i)
    {
        dr = chain.pos[i] + chain.pos[i + 1];
        inv_r = dr.reNorm();
        inv_r3 = (inv_r * inv_r * inv_r);
        velIndepAcc[chainIndex[i]]     += dr * (inv_r3 * this->m[chainIndex[i + 2]]);
        velIndepAcc[chainIndex[i + 2]] -= dr * (inv_r3 * this->m[chainIndex[i]]);
    }

    for(size_t i = 0 ; i < size() ; ++i)
    {
        for(size_t j = i + 3 ; j < size(); ++j)
        {
            dr = this->dynState.pos[chainIndex[j]] - this->dynState.pos[chainIndex[i]];
            inv_r = dr.reNorm();
            inv_r3 = (inv_r * inv_r * inv_r);
            velIndepAcc[chainIndex[i]] += dr * (inv_r3 * this->m[chainIndex[j]]);
            velIndepAcc[chainIndex[j]] -= dr * (inv_r3 * this->m[chainIndex[i]]);
        }
    }

    for(size_t i = 0 ; i < size() - 1; ++i)
        this->acc[i] = velIndepAcc[chainIndex[i + 1]] - velIndepAcc[chainIndex[i]];
}
#endif

