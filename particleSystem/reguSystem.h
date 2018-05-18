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
template<typename Interaction, typename EvolvedData, typename Regularitor>
class reguSystem : public particleSystem<reguSystem<Interaction, EvolvedData, Regularitor>, EvolvedData>
{
public:
    /////////////////////////////////Type Define////////////////////////////////////
    typedef typename EvolvedData::Scalar              Scalar;
    typedef typename EvolvedData::Vector              Vector;
    typedef typename EvolvedData::VectorArray         VectorArray;
    typedef typename EvolvedData::ScalarArray         ScalarArray;
    typedef std::array<Scalar, EvolvedData::volume()> PlainArray;
    /////////////////////////////////Interface /////////////////////////////////////
    constexpr static size_t size(){return EvolvedData::size();}
    constexpr static size_t volume(){return EvolvedData::volume();}
    void  advancePos(Scalar timeStepSize);
    void  advanceVel(Scalar timeStepSize);
    const reguSystem& operator=(const reguSystem& other);
    //////////////////////Base class Interface implement////////////////////////////
    std::istream& read(std::istream&);
    void   load(PlainArray& data);
    Scalar timeScale(Scalar scale);
    ///////////////////////////////Member variables/////////////////////////////////
private:
#ifdef KAHAN_SUMMATION
    EvolvedData roundoffErr;
#endif
    VectorArray velIndepAcc;
    VectorArray velDepAcc;
    Interaction velDepForce;
    Regularitor regular;
    /////////////////////////////private function///////////////////////////////////
private:
      void advanceOmega(Scalar stepSize);
      void advanceB(Scalar stepSize);
      void kickVel(Scalar stepSize);
      void kickAuxiVel(Scalar stepSize);
      void updateAccWith(VectorArray& vel);
      void updateVelIndepAcc();
};
////////////////////////////implement function//////////////////////////////////
template<typename Interaction, typename EvolvedData, typename Regularitor>
const reguSystem<Interaction, EvolvedData, Regularitor>& reguSystem<Interaction, EvolvedData, Regularitor>::operator=(const reguSystem& other)
{
    this->dynState    = other.dynState;
    this->acc         = other.acc;
    this->mass        = other.mass;
    this->radius      = other.radius;
    this->type        = other.type;
    this->velIndepAcc = other.velIndepAcc;
    this->velDepAcc   = other.velDepAcc;
    return *this;
}

template<typename Interaction, typename EvolvedData, typename Regularitor>
void reguSystem<Interaction, EvolvedData, Regularitor>::advancePos(Scalar timeStepSize)
{
    Scalar physicalTime = regular.getPhysicalPosTime(this->mass, this->dynState, timeStepSize);
#ifdef KAHAN_SUMMATION
    KahanAdvance(this->pos, this->vel, roundoffErr.pos, physicalTime);
#else
    advanceVariable(this->pos, this->vel, physicalTime);
#endif
    this->time += physicalTime;
}

template<typename Interaction, typename EvolvedData, typename Regularitor>
void reguSystem<Interaction, EvolvedData, Regularitor>::advanceVel(Scalar timeStepSize)
{
    Scalar physicalTime = regular.getPhysicalVelTime(this->mass, this->dynState, timeStepSize);
    Scalar halfTime = 0.5*physicalTime;
    updateVelIndepAcc();
    
    updateAccWith(this->vel);
    kickAuxiVel(halfTime);
    
    updateAccWith(this->dynState.auxiVel);
    advanceB(physicalTime);
    advanceOmega(physicalTime);
    kickVel(physicalTime);
    
    updateAccWith(this->vel);
    kickAuxiVel(halfTime);
}

template<typename Interaction, typename EvolvedData, typename Regularitor>
std::istream& reguSystem<Interaction, EvolvedData, Regularitor>::read(std::istream& input)
{
    input >> this->time;
    size_t id;
    for(size_t i = 0 ; i < size() ; ++i)
        input >> id >> this->type[i] >> this->mass[i] >> this->radius[i] >> this->pos[i] >> this->vel[i];
    memset(&(this->acc[0]), 0, sizeof(Vector)*size());
    MoveToCentralMassCoordinate(this->mass, this->pos);
    MoveToCentralMassCoordinate(this->mass, this->vel);
    this->dynState.initAddiVariable(this->mass);
#ifdef KAHAN_SUMMATION
    roundoffErr.setZero();
#endif
    return input;
}

template<typename Interaction, typename EvolvedData, typename Regularitor>
void reguSystem<Interaction, EvolvedData, Regularitor>::load(PlainArray& data)
{
    this->dynState.array() = data;
    this->dynState.auxiVel = this->vel;
#ifdef KAHAN_SUMMATION
    roundoffErr.setZero();
#endif
}

template<typename Interaction, typename EvolvedData, typename Regularitor>
typename EvolvedData::Scalar reguSystem<Interaction, EvolvedData, Regularitor>::timeScale(Scalar scale)
{
    return regular.getPhysicalPosTime(this->mass, this->dynState, scale);
}
/////////////////////////////private function///////////////////////////////////
template<typename Interaction, typename EvolvedData, typename Regularitor>
void reguSystem<Interaction, EvolvedData, Regularitor>::kickVel(Scalar stepSize)
{
#ifdef KAHAN_SUMMATION
    KahanAdvance(this->vel, this->acc, roundoffErr.vel, stepSize);
#else
    advanceVariable(this->vel, this->acc, stepSize);
#endif
}

template<typename Interaction, typename EvolvedData, typename Regularitor>
void reguSystem<Interaction, EvolvedData, Regularitor>::kickAuxiVel(Scalar stepSize)
{
#ifdef KAHAN_SUMMATION
    KahanAdvance(this->dynState.auxiVel, this->acc, roundoffErr.auxiVel, stepSize);
#else
    advanceVariable(this->dynState.auxiVel, this->acc, stepSize);
#endif
}

template<typename Interaction, typename EvolvedData, typename Regularitor>
void reguSystem<Interaction, EvolvedData, Regularitor>::advanceOmega(Scalar stepSize)
{
    Scalar dOmega = 0;
    for(size_t i = 0 ; i < size() ; ++i)
        dOmega += (this->velIndepAcc[i]*this->dynState.auxiVel[i])*(this->mass[i]);
    
#ifdef KAHAN_SUMMATION
    KahanAdvance(this->dynState.omega, dOmega*stepSize, roundoffErr.omega);
#else
    this->dynState.omega += dOmega*stepSize;
#endif
}

template<typename Interaction, typename EvolvedData, typename Regularitor>
void reguSystem<Interaction, EvolvedData, Regularitor>::advanceB(Scalar stepSize)
{
    Scalar dE = 0;
    for(size_t i = 0 ; i < size() ; ++i)
        dE -= (this->velDepAcc[i]*this->dynState.auxiVel[i])*(this->mass[i]);
    
#ifdef KAHAN_SUMMATION
    KahanAdvance(this->dynState.bindE, dE*stepSize, roundoffErr.bindE);
#else
    this->dynState.bindE += dE*stepSize;
#endif
}

template<typename Interaction, typename EvolvedData, typename Regularitor>
void reguSystem<Interaction, EvolvedData, Regularitor>::updateVelIndepAcc()
{
    Vector dr(0.0, 0.0, 0.0);
    Scalar inv_r  = 1;
    Scalar inv_r3 = 1;
    memset(&(velIndepAcc[0]), 0, sizeof(Vector)*size());
    for(size_t i = 0 ; i < size() ; ++i)
    {
        for(size_t j = i + 1 ; j < size() ; ++j)
        {
            dr     = this->pos[j] - this->pos[i];
            inv_r  = dr.reNorm();
            inv_r3 = inv_r*inv_r*inv_r;
            velIndepAcc[i] += dr*(inv_r3*this->mass[j]);
            velIndepAcc[j] -= dr*(inv_r3*this->mass[i]);
        }
    }
}

template<typename Interaction, typename EvolvedData, typename Regularitor>
void reguSystem<Interaction, EvolvedData, Regularitor>::updateAccWith(VectorArray& velocity)
{
    memset(&(velDepAcc[0]), 0, sizeof(Vector)*size());
    Vector dr;
    Vector dv;
    for(size_t i = 0 ; i < size() ; ++i)
    {
        for(size_t j = i + 1 ; j < size() ; ++j)
        {
            dr = this->pos[i] - this->pos[j];
            dv = velocity[i] - velocity[j];
            velDepForce(this->mass[i], this->mass[j], dr, dv,
                        velocity[i], velocity[j],
                        velDepAcc[i], velDepAcc[j]);
        }
    }
    for(size_t i = 0 ; i < size() ; ++i)
        this->acc[i] = velIndepAcc[i] + velDepAcc[i];
}

/*==============================Specialization================================*/

template<typename EvolvedData, typename Regularitor>
class reguSystem<Newtonian<typename EvolvedData::Scalar>, EvolvedData, Regularitor> :
public particleSystem<reguSystem<Newtonian<typename EvolvedData::Scalar>, EvolvedData, Regularitor>, EvolvedData>
{
public:
    /////////////////////////////////size_terface /////////////////////////////////////
    typedef typename EvolvedData::Scalar        Scalar;
    typedef typename EvolvedData::Vector        Vector;
    typedef typename EvolvedData::VectorArray   VectorArray;
    typedef typename EvolvedData::ScalarArray   ScalarArray;
    typedef std::array<Scalar, EvolvedData::volume()> PlainArray;
    
    constexpr static size_t size(){return EvolvedData::size();}
    constexpr static size_t volume(){return EvolvedData::volume();}
    void advancePos(Scalar timeStepSize);
    void advanceVel(Scalar timeStepSize);
    //////////////////////Base class size_terface implement////////////////////////////
    std::istream& read(std::istream&);
    void load(PlainArray& data);
    Scalar timeScale(Scalar scale);
    ///////////////////////////////Member variables/////////////////////////////////
private:
#ifdef KAHAN_SUMMATION
    EvolvedData  roundoffErr;
#endif
    Regularitor regular;
    /////////////////////////////private function///////////////////////////////////
private:
    void advanceOmega(Scalar stepSize);
    void updateVelIndepAcc();
};
////////////////////////////implement function//////////////////////////////////

template<typename EvolvedData, typename Regularitor>
void reguSystem<Newtonian<typename EvolvedData::Scalar>, EvolvedData, Regularitor>::advancePos(Scalar timeStepSize)
{
    Scalar physicalTime = regular.getPhysicalPosTime(this->mass, this->dynState, timeStepSize);
#ifdef KAHAN_SUMMATION
    KahanAdvance(this->pos, this->vel, roundoffErr.pos, physicalTime);
#else
    advanceVariable(this->pos, this->vel, physicalTime);
#endif
    this->time += physicalTime;
}

template<typename EvolvedData, typename Regularitor>
void reguSystem<Newtonian<typename EvolvedData::Scalar>, EvolvedData, Regularitor>::advanceVel(Scalar timeStepSize)
{
    Scalar physicalTime = regular.getPhysicalVelTime(this->mass, this->dynState, timeStepSize);
    updateVelIndepAcc();
    advanceOmega(0.5*physicalTime);
#ifdef KAHAN_SUMMATION
    KahanAdvance(this->vel, this->acc, roundoffErr.vel, physicalTime);
#else
    advanceVariable(this->vel, this->acc, physicalTime);
#endif
    advanceOmega(0.5*physicalTime);
}

template<typename EvolvedData, typename Regularitor>
std::istream& reguSystem<Newtonian<typename EvolvedData::Scalar>, EvolvedData, Regularitor>::read(std::istream& input)
{
    input >> this->time;
    size_t id;
    for(size_t i = 0 ; i < size() ; ++i)
        input >> id >> this->type[i] >> this->mass[i] >> this->radius[i] >> this->pos[i] >> this->vel[i];
    memset(&(this->acc[0]), 0, sizeof(Vector)*size());
    MoveToCentralMassCoordinate(this->mass, this->pos);
    MoveToCentralMassCoordinate(this->mass, this->vel);
    this->dynState.initAddiVariable(this->mass);
#ifdef KAHAN_SUMMATION
    roundoffErr.setZero();
#endif
    return input;
}

template<typename EvolvedData, typename Regularitor>
void reguSystem<Newtonian<typename EvolvedData::Scalar>, EvolvedData, Regularitor>::load(PlainArray& data)
{
    this->dynState.array() = data;
#ifdef KAHAN_SUMMATION
    roundoffErr.setZero();
#endif
}

template<typename EvolvedData, typename Regularitor>
typename EvolvedData::Scalar reguSystem<Newtonian<typename EvolvedData::Scalar>, EvolvedData, Regularitor>::timeScale(Scalar scale)
{
    return regular.getPhysicalPosTime(this->mass, this->dynState, scale);
}
/////////////////////////////private function///////////////////////////////////
template<typename EvolvedData, typename Regularitor>
void reguSystem<Newtonian<typename EvolvedData::Scalar>, EvolvedData, Regularitor>::advanceOmega(Scalar stepSize)
{
    Scalar dOmega = 0;
    for(size_t i = 0 ; i < size() ; ++i)
        dOmega += (this->acc[i]*this->vel[i])*(this->mass[i]);

#ifdef KAHAN_SUMMATION
    KahanAdvance(this->dynState.omega, dOmega*stepSize, roundoffErr.omega);
#else
    this->dynState.omega += dOmega*stepSize;
#endif
}

template<typename EvolvedData, typename Regularitor>
void reguSystem<Newtonian<typename EvolvedData::Scalar>, EvolvedData, Regularitor>::updateVelIndepAcc()
{
    Vector dr(0.0, 0.0, 0.0);
    Scalar inv_r  = 1;
    Scalar inv_r3 = 1;
    std::for_each(this->acc.begin(), this->acc.end(), [](vec3<Scalar>& v){ v.setZero(); });
    for(size_t i = 0 ; i < size() ; ++i)
    {
        for(size_t j = i + 1 ; j < size() ; ++j)
        {
            dr     = this->pos[j] - this->pos[i];
            inv_r  = dr.reNorm();
            inv_r3 = inv_r*inv_r*inv_r;
            this->acc[i] += dr*(inv_r3*this->mass[j]);
            this->acc[j] -= dr*(inv_r3*this->mass[i]);
        }
    }
}
#endif

