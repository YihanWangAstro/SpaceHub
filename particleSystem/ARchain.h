////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Filename:chainSystem.h                                                                                              //
//Author:Yihan Wang                                                                                                   //
//                                                                                                                    //
//                                                                                                                    //
//Description:                                                                                                        //
//  The main class of this n-body code. This class includes                                                           //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef ARCHAIN_H
#define ARCHAIN_H
#include "../particleSystem.h"
#include "../libs.h"
#include "chain.h"
#include <cstring>
#include <algorithm>
template<typename Interaction, typename EvolvedData, typename Regularitor>
class ARchain : public particleSystem<ARchain<Interaction, EvolvedData, Regularitor>, EvolvedData>
{
public:
    /////////////////////////////////Type Define////////////////////////////////////
    typedef typename EvolvedData::Scalar              Scalar;
    typedef typename EvolvedData::Vector              Vector;
    typedef typename EvolvedData::VectorArray         VectorArray;
    typedef typename EvolvedData::ScalarArray         ScalarArray;
    typedef std::array<size_t, EvolvedData::size()>    IndexArray;
    typedef std::array<Scalar, EvolvedData::volume()> PlainArray;
    /////////////////////////////////Interface /////////////////////////////////////
    constexpr static size_t size()
    {
        return EvolvedData::size();
    }
    void advancePos(Scalar timeStepSize);
    void advanceVel(Scalar timeStepSize);
    const ARchain& operator=(const ARchain& other);
    //////////////////////Base class Interface implement////////////////////////////
    std::istream& read(std::istream&);
    void load(PlainArray& data);
    Scalar timeScale(Scalar scale);
    PlainArray&   array()
    {
        return chain.array();
    }
    ///////////////////////////////Member variables/////////////////////////////////
private:
#ifdef KAHAN_SUMMATION
    EvolvedData roundoffErr;
#endif
    EvolvedData chain;
    IndexArray  chainIndex;
    VectorArray velIndepAcc;
    VectorArray chainVelIndepAcc;
    VectorArray velDepAcc;
    VectorArray chainVelDepAcc;
    Interaction velDepForce;
    Regularitor regular;
    /////////////////////////////private function///////////////////////////////////
private:
    void advanceOmega(Scalar stepSize);
    void advanceB(Scalar stepSize);
    void kickVel(Scalar stepSize);
    void kickAuxiVel(Scalar stepSize);
    void updateAccWith(VectorArray& vel, VectorArray& chainVel);
    void updateVelIndepAcc();
    void updateChain();
};
////////////////////////////implement function//////////////////////////////////
template<typename Interaction, typename EvolvedData, typename Regularitor>
const ARchain<Interaction, EvolvedData, Regularitor>& ARchain<Interaction, EvolvedData, Regularitor>::operator=
(const ARchain& other)
{
    this->dynState         = other.dynState;
    this->acc              = other.acc;
    this->mass             = other.mass;
    this->radius           = other.radius;
    this->type             = other.type;
    this->chain            = other.chain;
    this->chainIndex       = other.chainIndex;
    this->velIndepAcc      = other.velIndepAcc;
    this->chainVelIndepAcc = other.chainVelIndepAcc;
    this->velDepAcc        = other.velDepAcc;
    this->chainVelDepAcc   = other.chainVelDepAcc;
    return *this;
}

template<typename Interaction, typename EvolvedData, typename Regularitor>
void ARchain<Interaction, EvolvedData, Regularitor>::advancePos(Scalar timeStepSize)
{
    Scalar physicalTime = regular.getPhysicalPosTime(this->mass, this->dynState, timeStepSize);
#ifdef KAHAN_SUMMATION
    KahanAdvance(chain.pos, chain.vel, roundoffErr.pos, physicalTime);
#else
    advanceVariable(chain.pos, chain.vel, physicalTime);
#endif
    chain::synCartesian(chain.pos, this->pos, chainIndex);
    MoveToCentralMassCoordinate(this->mass, this->pos);
    this->chain.time += physicalTime;
    this->dynState.time = this->chain.time;
}

template<typename Interaction, typename EvolvedData, typename Regularitor>
void ARchain<Interaction, EvolvedData, Regularitor>::advanceVel(Scalar timeStepSize)
{
    Scalar physicalTime = regular.getPhysicalVelTime(this->mass, this->dynState, timeStepSize);
    Scalar halfTime     = 0.5 * physicalTime;
    updateVelIndepAcc();
    updateAccWith(this->vel, chain.vel);
    kickAuxiVel(halfTime);
    updateAccWith(this->dynState.auxiVel, chain.auxiVel);
    advanceB(physicalTime);
    advanceOmega(physicalTime);
    kickVel(physicalTime);
    updateAccWith(this->vel, chain.vel);
    kickAuxiVel(halfTime);
}

template<typename Interaction, typename EvolvedData, typename Regularitor>
std::istream& ARchain<Interaction, EvolvedData, Regularitor>::read(std::istream& input)
{
    input >> this->time;
    size_t id;
    
    for(size_t i = 0 ; i < size() ; ++i)
        input >> id >> this->type[i] >> this->mass[i] >> this->radius[i] >> this->pos[i] >> this->vel[i];
        
    memset(&(this->acc[0]), 0, sizeof(Vector)*size());
    this->dynState.moveToCentralMassCoords(this->mass);
    this->dynState.initAddiVariable(this->mass);
    chain::getChainIndex(this->pos, chainIndex);
    this->dynState.toChain(chain, chainIndex);
#ifdef KAHAN_SUMMATION
    roundoffErr.setZero();
#endif
    return input;
}

template<typename Interaction, typename EvolvedData, typename Regularitor>
void ARchain<Interaction, EvolvedData, Regularitor>::load(PlainArray& data)
{
    chain.array() = data;
    // chain.auxiVel = chain.vel;
    chain.toCartesian(this->dynState, chainIndex);
    this->dynState.moveToCentralMassCoords(this->mass);
    updateChain();
#ifdef KAHAN_SUMMATION
    roundoffErr.setZero();
#endif
}

template<typename Interaction, typename EvolvedData, typename Regularitor>
typename EvolvedData::Scalar ARchain<Interaction, EvolvedData, Regularitor>::timeScale(Scalar scale)
{
    return regular.getPhysicalPosTime(this->mass, this->dynState, scale);
}
/////////////////////////////private function///////////////////////////////////
template<typename Interaction, typename EvolvedData, typename Regularitor>
void ARchain<Interaction, EvolvedData, Regularitor>::kickVel(Scalar stepSize)
{
#ifdef KAHAN_SUMMATION
    KahanAdvance(chain.vel, this->acc, roundoffErr.vel, stepSize);
#else
    advanceVariable(chain.vel, this->acc, stepSize);
#endif
    chain::synCartesian(chain.vel, this->vel, chainIndex);
    MoveToCentralMassCoordinate(this->mass, this->vel);
}

template<typename Interaction, typename EvolvedData, typename Regularitor>
void ARchain<Interaction, EvolvedData, Regularitor>::kickAuxiVel(Scalar stepSize)
{
#ifdef KAHAN_SUMMATION
    KahanAdvance(chain.auxiVel, this->acc, roundoffErr.auxiVel, stepSize);
#else
    advanceVariable(chain.auxiVel, this->acc, stepSize);
#endif
    chain::synCartesian(chain.auxiVel, this->dynState.auxiVel, chainIndex);
    MoveToCentralMassCoordinate(this->mass, this->dynState.auxiVel);
}

template<typename Interaction, typename EvolvedData, typename Regularitor>
void ARchain<Interaction, EvolvedData, Regularitor>::advanceOmega(Scalar stepSize)
{
    Scalar dOmega = 0;
    
    for(size_t i = 0 ; i < size() ; ++i)
        dOmega += (this->velIndepAcc[i] * this->dynState.auxiVel[i]) * (this->mass[i]);
        
#ifdef KAHAN_SUMMATION
    KahanAdvance(this->chain.omega, dOmega * stepSize, roundoffErr.omega);
#else
    this->chain.omega += dOmega * stepSize;
#endif
    this->dynState.omega = this->chain.omega;
}

template<typename Interaction, typename EvolvedData, typename Regularitor>
void ARchain<Interaction, EvolvedData, Regularitor>::advanceB(Scalar stepSize)
{
    Scalar dE = 0;
    
    for(size_t i = 0 ; i < size() ; ++i)
        dE -= (this->velDepAcc[i] * this->dynState.auxiVel[i]) * (this->mass[i]);
        
#ifdef KAHAN_SUMMATION
    KahanAdvance(this->chain.bindE, dE * stepSize, roundoffErr.bindE);
#else
    this->chain.bindE += dE * stepSize;
#endif
    this->dynState.bindE = this->chain.bindE;
}

template<typename Interaction, typename EvolvedData, typename Regularitor>
void ARchain<Interaction, EvolvedData, Regularitor>::updateChain()
{
    IndexArray newIndex;
    chain::getChainIndex(this->pos, newIndex);
    
    if(chain::IsDiff(chainIndex, newIndex))
    {
        chain::updateChain(this->chain.pos, chainIndex,
                           newIndex);//Update chain index and update new position chain with old chain
        chainIndex = newIndex;
        chain::synChain(this->vel, chain.vel, chainIndex);
        chain::synChain(this->dynState.auxiVel, chain.auxiVel, chainIndex);
    }
}

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
        velIndepAcc[chainIndex[i]]     += dr * (inv_r3 * this->mass[chainIndex[i + 1]]);
        velIndepAcc[chainIndex[i + 1]] -= dr * (inv_r3 * this->mass[chainIndex[i]]);
    }
    
    for(size_t i = 0 ; i < size() - 2; ++i)
    {
        dr = chain.pos[i] + chain.pos[i + 1];
        inv_r = dr.reNorm();
        inv_r3 = (inv_r * inv_r * inv_r);
        velIndepAcc[chainIndex[i]]     += dr * (inv_r3 * this->mass[chainIndex[i + 2]]);
        velIndepAcc[chainIndex[i + 2]] -= dr * (inv_r3 * this->mass[chainIndex[i]]);
    }
    
    for(size_t i = 0 ; i < size() ; ++i)
    {
        for(size_t j = i + 3 ; j < size(); ++j)
        {
            dr = this->pos[chainIndex[j]] - this->pos[chainIndex[i]];
            inv_r = dr.reNorm();
            inv_r3 = (inv_r * inv_r * inv_r);
            velIndepAcc[chainIndex[i]] += dr * (inv_r3 * this->mass[chainIndex[j]]);
            velIndepAcc[chainIndex[j]] -= dr * (inv_r3 * this->mass[chainIndex[i]]);
        }
    }
    
    for(size_t i = 0 ; i < size() - 1; ++i)
        chainVelIndepAcc[i] = velIndepAcc[chainIndex[i + 1]] - velIndepAcc[chainIndex[i]];
        
    /*for(size_t i = 0 ; i < size() - 1; ++i)
        this->acc[i] = velIndepAcc[chainIndex[i + 1]] - velIndepAcc[chainIndex[i]];*/
}

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
        velDepForce(this->mass[chainIndex[i]], this->mass[chainIndex[i - 1]], dr, dv,
                    velocity[chainIndex[i]], velocity[chainIndex[i - 1]],
                    velDepAcc[chainIndex[i]], velDepAcc[chainIndex[i - 1]]);
    }
    
    for(int i = 2 ; i < size(); ++i)//update chain data(second near pair)
    {
        dr = chain.pos[i - 2] + chain.pos[i - 1];
        dv = chainVel[i - 2] + chainVel[i - 1];
        velDepForce(this->mass[chainIndex[i]], this->mass[chainIndex[i - 2]], dr, dv,
                    velocity[chainIndex[i]], velocity[chainIndex[i - 2]],
                    velDepAcc[chainIndex[i]], velDepAcc[chainIndex[i - 2]]);
    }
    
    for(int i = 0 ; i < size(); ++i)//update non-chain data
    {
        for(int j = 0 ; j < i - 2; ++j)
        {
            dr = this->pos[chainIndex[i]] - this->pos[chainIndex[j]];
            dv = velocity[chainIndex[i]] - velocity[chainIndex[j]];
            velDepForce(this->mass[chainIndex[i]], this->mass[chainIndex[j]], dr, dv,
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

template<typename EvolvedData, typename Regularitor>
class ARchain<Newtonian<typename EvolvedData::Scalar>, EvolvedData, Regularitor> :
    public particleSystem<ARchain<Newtonian<typename EvolvedData::Scalar>, EvolvedData, Regularitor>, EvolvedData>
{
public:
    /////////////////////////////////Type Define////////////////////////////////////
    typedef typename EvolvedData::Scalar              Scalar;
    typedef typename EvolvedData::Vector              Vector;
    typedef typename EvolvedData::VectorArray         VectorArray;
    typedef typename EvolvedData::ScalarArray         ScalarArray;
    typedef std::array<size_t, EvolvedData::size()>    IndexArray;
    typedef std::array<Scalar, EvolvedData::volume()> PlainArray;
    /////////////////////////////////Interface /////////////////////////////////////
    constexpr static size_t size()
    {
        return EvolvedData::size();
    }
    void advancePos(Scalar timeStepSize);
    void advanceVel(Scalar timeStepSize);
    const ARchain& operator=(const ARchain& other);
    //////////////////////Base class Interface implement////////////////////////////
    std::istream& read(std::istream&);
    void load(PlainArray& data);
    Scalar timeScale(Scalar scale);
    PlainArray&   array()
    {
        return chain.array();
    }
    ///////////////////////////////Member variables/////////////////////////////////
private:
#ifdef KAHAN_SUMMATION
    EvolvedData roundoffErr;
#endif
    EvolvedData chain;
    IndexArray  chainIndex;
    VectorArray velIndepAcc;
    Regularitor regular;
    /////////////////////////////private function///////////////////////////////////
private:
    void advanceOmega(Scalar stepSize);
    void updateVelIndepAcc();
    void updateChain();
};
////////////////////////////implement function//////////////////////////////////
template<typename EvolvedData, typename Regularitor>
const ARchain<Newtonian<typename EvolvedData::Scalar>, EvolvedData, Regularitor>&
ARchain<Newtonian<typename EvolvedData::Scalar>, EvolvedData, Regularitor>::operator=(const ARchain& other)
{
    this->dynState         = other.dynState;
    this->acc              = other.acc;
    this->mass             = other.mass;
    this->radius           = other.radius;
    this->type             = other.type;
    this->chain            = other.chain;
    this->chainIndex       = other.chainIndex;
    this->velIndepAcc      = other.velIndepAcc;
    return *this;
}

template<typename EvolvedData, typename Regularitor>
void ARchain<Newtonian<typename EvolvedData::Scalar>, EvolvedData, Regularitor>::advancePos(Scalar timeStepSize)
{
    Scalar physicalTime = regular.getPhysicalPosTime(this->mass, this->dynState, timeStepSize);
#ifdef KAHAN_SUMMATION
    KahanAdvance(chain.pos, chain.vel, roundoffErr.pos, physicalTime);
#else
    advanceVariable(chain.pos, chain.vel, physicalTime);
#endif
    chain::synCartesian(chain.pos, this->pos, chainIndex);
    MoveToCentralMassCoordinate(this->mass, this->pos);
    this->chain.time += physicalTime;
    this->dynState.time = this->chain.time;
}

template<typename EvolvedData, typename Regularitor>
void ARchain<Newtonian<typename EvolvedData::Scalar>, EvolvedData, Regularitor>::advanceVel(Scalar timeStepSize)
{
    Scalar physicalTime = regular.getPhysicalVelTime(this->mass, this->dynState, timeStepSize);
    Scalar halfTime     = 0.5 * physicalTime;
    updateVelIndepAcc();
    advanceOmega(halfTime);
#ifdef KAHAN_SUMMATION
    KahanAdvance(chain.vel, this->acc, roundoffErr.vel, physicalTime);
#else
    advanceVariable(chain.vel, this->acc, physicalTime);
#endif
    chain::synCartesian(chain.vel, this->vel, chainIndex);
    MoveToCentralMassCoordinate(this->mass, this->vel);
    advanceOmega(halfTime);
}

template<typename EvolvedData, typename Regularitor>
std::istream& ARchain<Newtonian<typename EvolvedData::Scalar>, EvolvedData, Regularitor>::read(std::istream& input)
{
    input >> this->time;
    size_t id;
    
    for(size_t i = 0 ; i < size() ; ++i)
        input >> id >> this->type[i] >> this->mass[i] >> this->radius[i] >> this->pos[i] >> this->vel[i];
        
    memset(&(this->acc[0]), 0, sizeof(Vector)*size());
    this->dynState.moveToCentralMassCoords(this->mass);
    this->dynState.initAddiVariable(this->mass);
    chain::getChainIndex(this->pos, chainIndex);
    this->dynState.toChain(chain, chainIndex);
#ifdef KAHAN_SUMMATION
    roundoffErr.setZero();
#endif
    return input;
}

template<typename EvolvedData, typename Regularitor>
void ARchain<Newtonian<typename EvolvedData::Scalar>, EvolvedData, Regularitor>::load(PlainArray& data)
{
    chain.array() = data;
    chain.toCartesian(this->dynState, chainIndex);
    this->dynState.moveToCentralMassCoords(this->mass);
    updateChain();
#ifdef KAHAN_SUMMATION
    roundoffErr.setZero();
#endif
}

template<typename EvolvedData, typename Regularitor>
typename EvolvedData::Scalar ARchain<Newtonian<typename EvolvedData::Scalar>, EvolvedData, Regularitor>::timeScale(
    Scalar scale)
{
    return regular.getPhysicalPosTime(this->mass, this->dynState, scale);
}
/////////////////////////////private function///////////////////////////////////
template<typename EvolvedData, typename Regularitor>
void ARchain<Newtonian<typename EvolvedData::Scalar>, EvolvedData, Regularitor>::advanceOmega(Scalar stepSize)
{
    Scalar dOmega = 0;
    
    for(size_t i = 0 ; i < size() ; ++i)
        dOmega += (this->velIndepAcc[i] * this->dynState.vel[i]) * (this->mass[i]);
        
#ifdef KAHAN_SUMMATION
    KahanAdvance(this->chain.omega, dOmega * stepSize, roundoffErr.omega);
#else
    this->chain.omega += dOmega * stepSize;
#endif
    this->dynState.omega = this->chain.omega;
}

template<typename EvolvedData, typename Regularitor>
void ARchain<Newtonian<typename EvolvedData::Scalar>, EvolvedData, Regularitor>::updateChain()
{
    IndexArray newIndex;
    chain::getChainIndex(this->pos, newIndex);
    
    if(chain::IsDiff(chainIndex, newIndex))
    {
        chain::updateChain(this->chain.pos, chainIndex,
                           newIndex);//Update chain index and update new position chain with old chain
        chainIndex = newIndex;
        chain::synChain(this->vel, chain.vel, chainIndex);
    }
}

template<typename EvolvedData, typename Regularitor>
void ARchain<Newtonian<typename EvolvedData::Scalar>, EvolvedData, Regularitor>::updateVelIndepAcc()
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
        velIndepAcc[chainIndex[i]]     += dr * (inv_r3 * this->mass[chainIndex[i + 1]]);
        velIndepAcc[chainIndex[i + 1]] -= dr * (inv_r3 * this->mass[chainIndex[i]]);
    }
    
    for(size_t i = 0 ; i < size() - 2; ++i)
    {
        dr = chain.pos[i] + chain.pos[i + 1];
        inv_r = dr.reNorm();
        inv_r3 = (inv_r * inv_r * inv_r);
        velIndepAcc[chainIndex[i]]     += dr * (inv_r3 * this->mass[chainIndex[i + 2]]);
        velIndepAcc[chainIndex[i + 2]] -= dr * (inv_r3 * this->mass[chainIndex[i]]);
    }
    
    for(size_t i = 0 ; i < size() ; ++i)
    {
        for(size_t j = i + 3 ; j < size(); ++j)
        {
            dr = this->pos[chainIndex[j]] - this->pos[chainIndex[i]];
            inv_r = dr.reNorm();
            inv_r3 = (inv_r * inv_r * inv_r);
            velIndepAcc[chainIndex[i]] += dr * (inv_r3 * this->mass[chainIndex[j]]);
            velIndepAcc[chainIndex[j]] -= dr * (inv_r3 * this->mass[chainIndex[i]]);
        }
    }
    
    for(size_t i = 0 ; i < size() - 1; ++i)
        this->acc[i] = velIndepAcc[chainIndex[i + 1]] - velIndepAcc[chainIndex[i]];
}
#endif

