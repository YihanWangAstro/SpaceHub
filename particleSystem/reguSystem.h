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
template<typename Interaction, typename EvolvedData, typename Regularitor>
class reguSystem : public particleSystem<reguSystem<Interaction, EvolvedData, Regularitor>, EvolvedData>
{
public:

    typedef typename EvolvedData::Scalar              Scalar;
    typedef typename EvolvedData::Vector              Vector;
    typedef typename EvolvedData::VectorArray         VectorArray;
    typedef typename EvolvedData::ScalarArray         ScalarArray;
    typedef std::array<Scalar, EvolvedData::volume()> PlainArray;

    /** @brief Advance position one step with current velocity. */
    void  advancePos(Scalar timeStepSize);

    /** @brief Advance velocity one step with current acceleration. */
    void  advanceVel(Scalar timeStepSize);

    /** @brief Overload operator = */
    const reguSystem& operator=(const reguSystem& other);

    /** @brief Get the number of the particles.
     *  @return The particle number.
     *  @note  Overload base class size().
     */
    constexpr static size_t size()
    {
        return EvolvedData::size();
    }

    /** @brief Get the total dynamic scalar number.
     *  @return The dynamic scalar number.
     *  @note  Overload base class volume().
     */
    constexpr static size_t volume()
    {
        return EvolvedData::volume();
    }

    /** @brief Input data from standard c++ istream.
     *  @note  Overload base class read().
     */
    std::istream& read(std::istream&);

    /** @brief Load data from a plain array.
     *  @note  Overload base class load().
     */
    void   load(PlainArray& data);

    /** @brief Interface to rescale the time.
     *  @note  Overload base class timeScale().
     */
    Scalar timeScale(Scalar scale);

private:
#ifdef KAHAN_SUMMATION
    /**  @brief Co-evolved Data used by Kahan summation.
     *
     *   See details in https://en.wikipedia.org/wiki/Kahan_summation_algorithm .
     */
    EvolvedData roundoffErr;
#endif

    /** @brief Velocity independent acceleration array*/
    VectorArray velIndepAcc;

    /** @brief Velocity dependent acceleration array*/
    VectorArray velDepAcc;

    /** @brief Velocity dependent pair force functor*/
    Interaction velDepForce;

    /** @brief Regularization interface.*/
    Regularitor regular;

private:

    /** @brief Advance regularization variable omega.*/
    void advanceOmega(Scalar stepSize);

    /** @brief Advance regularization variable bindE.*/
    void advanceB(Scalar stepSize);

    /** @brief Advance velocity with current acceleration.*/
    void kickVel(Scalar stepSize);

    /** @brief Advance auxiliar velocity with current acceleration.*/
    void kickAuxiVel(Scalar stepSize);

    /** @brief Update velocity dependent acceleration with given velocity. Then update the total acceleration.*/
    void updateAccWith(VectorArray& vel);

    /** @brief Update velocity independent acceleration.*/
    void updateVelIndepAcc();
};

/** @brief Overload operator = */
template<typename Interaction, typename EvolvedData, typename Regularitor>
const reguSystem<Interaction, EvolvedData, Regularitor>& reguSystem<Interaction, EvolvedData, Regularitor>::operator=
(const reguSystem& other)
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

/**  @brief Advance position one step with current velocity.
  *
  *  Advance position array and physical time one step with current integration step size and velocity.
  *  @param  timeStepSize Integration step size, will be transfered to physical time in the function.
 */
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

/** @brief Advance velocity one step with current acceleration.
 *
 *  Advance velocity and auxiliary velocity array one step with current integration step size and accelerations.
 *  @param  timeStepSize Integration step size, will be transfered to physical time in the function.
 */
template<typename Interaction, typename EvolvedData, typename Regularitor>
void reguSystem<Interaction, EvolvedData, Regularitor>::advanceVel(Scalar timeStepSize)
{
    Scalar physicalTime = regular.getPhysicalVelTime(this->mass, this->dynState, timeStepSize);
    Scalar halfTime = 0.5 * physicalTime;
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

/** @brief Input data from standard c++ istream.
 *
 *  Implement of CRTP '>>' method.
 */
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

/** @brief Load data from a plain array.
 *
 *  Interface usded by integrator and ODE iterator. Load data from a plain array processed by itegrator and
 *  iterator. Derived class could overload this function to additional process.
 *
 *  @param data Plain scalar array.
 */
template<typename Interaction, typename EvolvedData, typename Regularitor>
void reguSystem<Interaction, EvolvedData, Regularitor>::load(PlainArray& data)
{
    this->dynState.array() = data;
    this->dynState.auxiVel = this->vel;
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
typename EvolvedData::Scalar reguSystem<Interaction, EvolvedData, Regularitor>::timeScale(Scalar scale)
{
    return regular.getPhysicalPosTime(this->mass, this->dynState, scale);
}

/** @brief Advance velocity with current acceleration.
 *
 *  Advance the velocity with current total acceleration variable 'acc'.
 *  @param stepSize Integration step size.
 */
template<typename Interaction, typename EvolvedData, typename Regularitor>
void reguSystem<Interaction, EvolvedData, Regularitor>::kickVel(Scalar stepSize)
{
#ifdef KAHAN_SUMMATION
    KahanAdvance(this->vel, this->acc, roundoffErr.vel, stepSize);
#else
    advanceVariable(this->vel, this->acc, stepSize);
#endif
}

/** @brief Advance auxiliar velocity with current acceleration.
 *
 *  Advance the auxiliar velocity with current total acceleration variable 'acc'.
 *  @param stepSize Integration step size.
 */
template<typename Interaction, typename EvolvedData, typename Regularitor>
void reguSystem<Interaction, EvolvedData, Regularitor>::kickAuxiVel(Scalar stepSize)
{
#ifdef KAHAN_SUMMATION
    KahanAdvance(this->dynState.auxiVel, this->acc, roundoffErr.auxiVel, stepSize);
#else
    advanceVariable(this->dynState.auxiVel, this->acc, stepSize);
#endif
}

/** @brief Advance regularization variable omega.
 *
 *  Advance omega with velocity independent acceleration and auxiliar velocity with physical time step size.
 *  @param stepSize Physical time step.
 */
template<typename Interaction, typename EvolvedData, typename Regularitor>
void reguSystem<Interaction, EvolvedData, Regularitor>::advanceOmega(Scalar stepSize)
{
    Scalar dOmega = 0;

    for(size_t i = 0 ; i < size() ; ++i)
        dOmega += (this->velIndepAcc[i] * this->dynState.auxiVel[i]) * (this->mass[i]);

#ifdef KAHAN_SUMMATION
    KahanAdvance(this->dynState.omega, dOmega * stepSize, roundoffErr.omega);
#else
    this->dynState.omega += dOmega * stepSize;
#endif
}

/** @brief Advance regularization variable bindE.
 *
 *  Advance bindE with velocity dependent acceleration and auxiliar velocity with physical time step size.
 *  @param stepSize Physical time step.
 */
template<typename Interaction, typename EvolvedData, typename Regularitor>
void reguSystem<Interaction, EvolvedData, Regularitor>::advanceB(Scalar stepSize)
{
    Scalar dE = 0;

    for(size_t i = 0 ; i < size() ; ++i)
        dE -= (this->velDepAcc[i] * this->dynState.auxiVel[i]) * (this->mass[i]);

#ifdef KAHAN_SUMMATION
    KahanAdvance(this->dynState.bindE, dE * stepSize, roundoffErr.bindE);
#else
    this->dynState.bindE += dE * stepSize;
#endif
}

/** @brief Update velocity independent acceleration.
 *
 *  Update velocity independent acceleration 'velIndepAcc' with Newtonian interaction.
 */
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
            inv_r3 = inv_r * inv_r * inv_r;
            velIndepAcc[i] += dr * (inv_r3 * this->mass[j]);
            velIndepAcc[j] -= dr * (inv_r3 * this->mass[i]);
        }
    }
}

/** @brief Update velocity dependent acceleration with given velocity. Then update the total acceleration.
 *
 *  Update the velocity dependent accelaration variable 'velDepAcc' with given velocity and velocity dependent pair
 *  force 'velDepForce'. Then update the total acceleration 'acc' by adding 'velIndepAcc' and 'velDepAcc'.
 */
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
/**  @brief Partial specialization of reguSystem for velocity independent system. */
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

    /** @brief Get the number of the particles.
     *  @return The particle number.
     *  @note  Overload base class size().
     */
    constexpr static size_t size()
    {
        return EvolvedData::size();
    }
    
    /** @brief Get the total dynamic scalar number.
     *  @return The dynamic scalar number.
     *  @note  Overload base class volume().
     */
    constexpr static size_t volume()
    {
        return EvolvedData::volume();
    }
    
    /** @brief Advance position one step with current velocity. */
    void advancePos(Scalar timeStepSize);
    
     /** @brief Advance velocity one step with current acceleration. */
    void advanceVel(Scalar timeStepSize);
    
    /** @brief Input data from standard c++ istream.
     *  @note  Overload base class read().
     */
    std::istream& read(std::istream&);
    
    /** @brief Load data from a plain array.
     *  @note  Overload base class load().
     */
    void load(PlainArray& data);
    
    /** @brief Interface to rescale the time.
     *  @note  Overload base class timeScale().
     */
    Scalar timeScale(Scalar scale);
    ///////////////////////////////Member variables/////////////////////////////////
private:
#ifdef KAHAN_SUMMATION
    /**  @brief Co-evolved Data used by Kahan summation.
     *
     *   See details in https://en.wikipedia.org/wiki/Kahan_summation_algorithm .
     */
    EvolvedData  roundoffErr;
#endif
    
    /** @brief Regularization interface.*/
    Regularitor regular;
    /////////////////////////////private function///////////////////////////////////
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

/** @brief Advance velocity one step with current acceleration.
 *
 *  Advance velocity array one step with current integration step size and accelerations.
 *  @param  timeStepSize Integration step size, will be transfered to physical time in the function.
 */
template<typename EvolvedData, typename Regularitor>
void reguSystem<Newtonian<typename EvolvedData::Scalar>, EvolvedData, Regularitor>::advanceVel(Scalar timeStepSize)
{
    Scalar physicalTime = regular.getPhysicalVelTime(this->mass, this->dynState, timeStepSize);
    updateVelIndepAcc();
    advanceOmega(0.5 * physicalTime);
#ifdef KAHAN_SUMMATION
    KahanAdvance(this->vel, this->acc, roundoffErr.vel, physicalTime);
#else
    advanceVariable(this->vel, this->acc, physicalTime);
#endif
    advanceOmega(0.5 * physicalTime);
}

/** @brief Input data from standard c++ istream.
 *
 *  Implement of CRTP '>>' method.
 */
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

/** @brief Load data from a plain array.
 *
 *  Interface usded by integrator and ODE iterator. Load data from a plain array processed by itegrator and
 *  iterator. Derived class could overload this function to additional process.
 *
 *  @param data Plain scalar array.
 */
template<typename EvolvedData, typename Regularitor>
void reguSystem<Newtonian<typename EvolvedData::Scalar>, EvolvedData, Regularitor>::load(PlainArray& data)
{
    this->dynState.array() = data;
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
typename EvolvedData::Scalar reguSystem<Newtonian<typename EvolvedData::Scalar>, EvolvedData, Regularitor>::timeScale(
    Scalar scale)
{
    return regular.getPhysicalPosTime(this->mass, this->dynState, scale);
}

/** @brief Advance regularization variable omega.
 *
 *  Advance omega with velocity independent acceleration and auxiliar velocity with physical time step size.
 *  @param stepSize Physical time step.
 */
template<typename EvolvedData, typename Regularitor>
void reguSystem<Newtonian<typename EvolvedData::Scalar>, EvolvedData, Regularitor>::advanceOmega(Scalar stepSize)
{
    Scalar dOmega = 0;

    for(size_t i = 0 ; i < size() ; ++i)
        dOmega += (this->acc[i] * this->vel[i]) * (this->mass[i]);

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
template<typename EvolvedData, typename Regularitor>
void reguSystem<Newtonian<typename EvolvedData::Scalar>, EvolvedData, Regularitor>::updateVelIndepAcc()
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
            dr     = this->pos[j] - this->pos[i];
            inv_r  = dr.reNorm();
            inv_r3 = inv_r * inv_r * inv_r;
            this->acc[i] += dr * (inv_r3 * this->mass[j]);
            this->acc[j] -= dr * (inv_r3 * this->mass[i]);
        }
    }
}
#endif

